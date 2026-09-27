#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import json
import sys
import os
import time
import subprocess
import mmap
import datetime
import platform
import faulthandler
import signal

from checkgms_utils import *
from checkgms_stable import *

#global states
first = True
if "Windows" in platform.system():
    rungms_path = "rungms.bat"
else:
    rungms_path = "misc/automation/rungms"
codecoverage = []

def run_test(input_file):
    global first
    global codecoverage

    filenum = input_file[0]
    input_file_path = input_file[1]
    input_file_path = re.sub("(\*|\?| |!|\$|#|&|\"|\'|\(|\)|\||<|>|\\\|;)", r"\\\1", input_file_path)
    short_input_file_path = input_file_path.split("/tests/", 1)[-1]

    rungms_path_tmp = run_arguments["path"]
    if rungms_path_tmp == "":
        rungms_path_tmp = rungms_path

    if run_arguments["bwrap"] and "Windows" not in platform.system():
        rungms_path_tmp = 'bwrap --setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}" --dev-bind / / --unshare-ipc ' + \
            rungms_path_tmp

    run_command = rungms_path_tmp + " " + input_file_path + " " + run_arguments["version"] + " " + \
        run_arguments["ncpus"]

    if "Windows" in platform.system():
        run_command = run_command + " > " + input_file_path.replace(".inp", ".log") + " 2>&1"
    else:
        run_command = run_command + " " + run_arguments["ncpus"] + \
        " > " + input_file_path.replace(".inp", ".log") + " 2>&1"

    if run_arguments["print_to_stdout"]:
        run_command = rungms_path_tmp + " " + input_file_path + " " + run_arguments["version"] + " " + \
            run_arguments["ncpus"]

        if "Windows" not in platform.system():
            run_command = run_command + " " + run_arguments["ncpus"]

    if run_arguments["pre"]:
        run_command = run_arguments["pre"] + " " + input_file_path + " && " + run_command

    if run_arguments["post"]:
        run_command = run_command + " && " + run_arguments["post"] + " " + input_file_path

    # TODO: Use Python 3's pathlib to juggle Windows / Unix paths
    if "Windows" in platform.system():
        run_command=run_command.replace("/","\\")

    if not run_arguments["coverage"]:
        if not run_arguments["no_counter"]:
            sys.stdout.write(
                c_box_small(
                    datetime.datetime.now().time()) + " " +
                l_box_small("Running input file") + " " +
                file_progress(
                    filenum,
                    len(input_file_paths)) + " " +
                short_input_file_path + "\n")

    if run_arguments["stderr"]:
        sys.stderr.write(
            l_box_small("Running input file") +
            " " +
            file_progress(
                filenum,
                len(input_file_paths)) +
            " " +
            short_input_file_path +
            "\n")
        sys.stderr.flush()

    if run_arguments["debug"] or run_arguments["dryrun"]:
        sys.stdout.write(run_command + "\n")

    if not run_arguments["dryrun"]:
        if "Windows" in platform.system():
            p=subprocess.Popen([run_command], shell=True, stdout=sys.stdout)
            try:
                p.wait(timeout=7200)
            except subprocess.TimeoutExpired:
                sys.stdout.write("Timeout expired for " + short_input_file_path + "\n")
                sys.stdout.flush()
                p.terminate()
                return
            except KeyboardInterrupt:
                sys.exit(1)
        else:
            p=subprocess.Popen([run_command], shell=True, stdout=sys.stdout, start_new_session=True)
            try:
                p.wait(timeout=7200)
            except subprocess.TimeoutExpired:
                sys.stdout.write("Timeout expired for " + short_input_file_path + "\n")
                sys.stdout.flush()
                os.killpg(os.getpgid(p.pid), signal.SIGTERM)
                return
            except KeyboardInterrupt:
                sys.exit(1)

    if run_arguments["coverage"]:
        try:
            result = subprocess.run(
                ["tests/coverage"],
                stdout=subprocess.PIPE).stdout.decode('utf-8')
        except KeyboardInterrupt:
            sys.exit(1)

        lines = result.splitlines()

        if first:
            sys.stdout.write(c_large("Test Input") + " " + c_small("% Total Coverage") + "\n")
            sys.stdout.write(seperator(columns=83) + "\n")
            codecoverage.append(["Source Files"])
            codecoverage.append(["Lines of Code"])

        codecoverage.append([short_input_file_path])

        for line in range(len(lines)):
            if "File" in lines[line]:
                sourcefile = lines[line].split("'")[-2]
                coverage = lines[line + 1].split(":")[-1].split("%")[0]
                linesofcode = lines[line + 1].split(" ")[-1]

                if "No executable lines" in lines[line + 1]:
                    coverage = "0.00"
                    linesofcode = 0

                if first:
                    codecoverage[0].append(sourcefile)
                    codecoverage[1].append(linesofcode)

                codecoverage[-1].append(coverage)

        if first:
            codecoverage[0].append("Total")
            loc_sum = 0
            for loc in codecoverage[1][1:-1]:
                loc_sum = loc_sum + int(loc)
            codecoverage[1].append(loc_sum)

        codecoverage[-1].append(lines[-1].split(":")[-1].split("%")[0])

        if first:
            first = False

        if float(codecoverage[-1][-1]) > float(codecoverage[-2][-1]):
            sys.stdout.write(l_box_large(
                codecoverage[-1][0]) + " " +
                r_small(codecoverage[-1][-1])+ "\n")
        else:
            if float(codecoverage[-2][-1]) < 100.0:
                sys.stdout.write(l_box_large(
                    codecoverage[-1][0]) + " " +
                    r_small(codecoverage[-1][-1]) + " " +
                    c_small("X") + "\n")
            else:
                sys.stdout.write(l_box_large(
                    codecoverage[-1][0]) + " " +
                    r_small(codecoverage[-1][-1]) + "\n")

if __name__ == '__main__':
    faulthandler.dump_traceback_later(12*3600, exit=True)
    input_files = []

    total_coverage = 0
    total_linesofcode = 0

    run_arguments = parse_arguments(checkgms=False)

    if run_arguments["supress"]== False:
        print(c_box("Run parameters"))
        print(json.dumps(run_arguments, indent=2))

    script_path = os.path.dirname(os.path.realpath(__file__))

    threads_num = 0
    try:
        threads_num = int(run_arguments["threads"])
    except ValueError:
        print("\"threads\" flag requires integer value, but " +
              run_arguments["threads"] + " is not an integer value!")
        sys.exit(1)

    if threads_num < 1:
        print("Why are you thinking that zero threads can do something?")
        sys.exit(1)

    if threads_num > 1:
        try:
            import concurrent.futures as cf
        except ImportError:
            print("Parallelization is not available with this python version!")
            sys.exit(1)

    if threads_num > 1 and run_arguments["coverage"]:
        print("coverage run can not parallelized!")
        sys.exit(1)

    # Populate the input_file_paths array
    input_file_paths = []
    input_file_paths = get_input_file_paths(
        filepath_string_match=run_arguments["filter_filepath"],
        folder_string_match=run_arguments["filter_folder"],
        file_string_match=run_arguments["filter_file"],
        folder_string_skip=run_arguments["skip_folder"],
        file_string_skip=run_arguments["skip_file"],
        script_path=script_path)

    # Loop through the input_file_paths array and validate
    for filenum, input_file_path in enumerate(input_file_paths, start=1):

        if ".inp" not in input_file_path[-4:]:
            continue

        short_input_file_path = input_file_path.split("/tests/", 1)[-1]

        # Check if user wants to skip input if log file already exists
        if run_arguments["skip_log"]:
            if (os.path.isfile(input_file_path.replace('.inp', '.log'))):
                print(
                    c_box_small(
                        datetime.datetime.now().time()),
                    l_box_small("Log file already exists!"),
                    file_progress(
                        filenum,
                        len(input_file_paths)),
                    short_input_file_path.replace('.inp', '.log'))
                continue

        # Search for prescence of "TRAVIS-CI SKIP"
        with open(input_file_path, 'r', encoding="utf-8", errors='ignore') as opened_file:
            parse_memory_map = mmap.mmap(
                opened_file.fileno(), 0, access=mmap.ACCESS_READ)

            if not run_arguments["no_skip"]:
                regex_string = "TRAVIS-CI SKIP"
                regex = re.compile(
                    str.encode(
                        regex_string,
                        'ascii'),
                    re.MULTILINE)
                match = regex.search(parse_memory_map)
                # If found then skip
                if match:
                    if run_arguments["debug"]:
                        print(
                            l_box_small("Skipping input file"),
                            file_progress(
                                filenum,
                                len(input_file_paths)),
                            short_input_file_path)
                    continue

            # If --test_type is passed in:
            if "small" in run_arguments["test_type"]:
                regex_string = "TRAVIS-CI SMALL"
            elif "medium" in run_arguments["test_type"]:
                regex_string = "TRAVIS-CI MEDIUM"
            elif "large" in run_arguments["test_type"]:
                regex_string = "TRAVIS-CI LARGE"
            elif "msucc" in run_arguments["test_type"]:
                regex_string = "TRAVIS-CI MSUCC"

            if len(run_arguments["test_type"]) > 0:
                regex = re.compile(
                    str.encode(
                        regex_string,
                        'ascii'),
                    re.MULTILINE)
                match = regex.search(parse_memory_map)
                if not match:
                    if run_arguments["debug"]:
                        print(
                            l_box_small("Skipping input file"),
                            file_progress(
                                filenum,
                                len(input_file_paths)),
                            short_input_file_path)
                    continue

        input_files.append([filenum, input_file_path])

    if threads_num == 1:
        for input_file in input_files:
            run_test(input_file)
    else:
        with cf.ThreadPoolExecutor(max_workers=threads_num) as executor:
            try:
                for _ in executor.map(run_test, input_files, chunksize=1):
                    pass
            except KeyboardInterrupt:
                sys.exit(1)

    if run_arguments["coverage"]:
        with open("tests/coverage.dat", "w") as coverage_file:
            for row in codecoverage:
                coverage_file.write(",".join([str(x) for x in row]))
                coverage_file.write("\n")
                pass

        print(seperator(columns=83))
        print(l_box_large("Total % Coverage"), r_small(codecoverage[-1][-1]))
        print(l_box_large("Total Lines of Code"), r_small(codecoverage[1][-1]))
