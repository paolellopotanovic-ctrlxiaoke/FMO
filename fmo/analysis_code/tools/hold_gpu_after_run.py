"""Reserve modest GPU memory immediately or after a production run succeeds."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def completion_state(run):
    records = {}
    for target in ("CDK8", "BACE"):
        for leg in ("solution", "bound"):
            path = run / "production" / target / "edge_000" / leg / "status.json"
            key = f"{target}/{leg}"
            try:
                records[key] = json.loads(path.read_text()).get("status", "UNKNOWN")
            except (FileNotFoundError, json.JSONDecodeError):
                records[key] = "PENDING"
    if all(value == "SUCCEEDED" for value in records.values()):
        return "READY", records
    if any(value.startswith(("FAILED", "INTERRUPTED")) for value in records.values()):
        return "FAILED", records
    return "WAITING", records


def write_status(directory, record):
    record["updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    record["watcher_pid"] = os.getpid()
    path = directory / "status.json"
    temporary = directory / "status.json.tmp"
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(path)


def reserve(memory_mib):
    import torch

    allocation = torch.empty(memory_mib * 1024 * 1024, dtype=torch.uint8, device="cuda:0")
    allocation.zero_()
    torch.cuda.synchronize()
    print(json.dumps({"status": "HOLDING", "pid": os.getpid(), "memory_mib": memory_mib,
                      "physical_device": os.environ.get("CUDA_VISIBLE_DEVICES"),
                      "continuous_compute": False}), flush=True)
    while allocation is not None:
        time.sleep(60)


def watch(run, directory, devices, memory_mib, interval, immediate=False):
    directory.mkdir(parents=True, exist_ok=True)
    children = []
    handles = []
    with (directory / "watcher.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            while not immediate:
                state, records = completion_state(run)
                write_status(directory, {"status": state, "legs": records,
                                        "devices": devices, "memory_mib_per_device": memory_mib})
                if state == "FAILED":
                    return 1
                if state == "READY":
                    break
                time.sleep(interval)
            if not immediate:
                time.sleep(10)
            for device in devices:
                environment = os.environ.copy()
                environment["CUDA_VISIBLE_DEVICES"] = str(device)
                handle = (directory / f"gpu_{device}.stdout").open("a")
                handles.append(handle)
                child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                    "--reserve", "--memory-mib", str(memory_mib)], env=environment,
                    stdout=handle, stderr=subprocess.STDOUT)
                children.append(child)
            write_status(directory, {"status": "STARTING_HOLDERS", "holder_pids": [child.pid for child in children],
                                    "devices": devices, "memory_mib_per_device": memory_mib})
            while True:
                if any(child.poll() is not None for child in children):
                    write_status(directory, {"status": "HOLDER_EXITED", "returncodes": [child.poll() for child in children]})
                    return 1
                initialized = all('"status": "HOLDING"' in (directory / f"gpu_{device}.stdout").read_text()
                                  for device in devices)
                if initialized:
                    write_status(directory, {"status": "HOLDING", "holder_pids": [child.pid for child in children],
                        "devices": devices, "memory_mib_per_device": memory_mib, "continuous_compute": False,
                        "start_policy": "immediate" if immediate else "after_production_success"})
                time.sleep(interval)
        finally:
            for child in children:
                if child.poll() is None:
                    child.terminate()
            for child in children:
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
            for handle in handles:
                handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--devices", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--interval", type=int, default=30)
    parser.add_argument("--reserve", action="store_true")
    parser.add_argument("--immediate", action="store_true")
    args = parser.parse_args()
    if args.memory_mib < 1 or args.interval < 1 or any(device < 0 for device in args.devices):
        parser.error("positive memory/interval and nonnegative devices required")
    if len(args.devices) != len(set(args.devices)):
        parser.error("duplicate devices are not allowed")
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
    if args.reserve:
        reserve(args.memory_mib)
    else:
        if args.run is None:
            parser.error("--run required for completion watcher")
        run = args.run.resolve(strict=True)
        raise SystemExit(watch(run, run / "gpu_hold", args.devices, args.memory_mib, args.interval,
                              immediate=args.immediate))


if __name__ == "__main__":
    main()
