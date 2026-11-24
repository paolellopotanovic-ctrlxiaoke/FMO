"""Build an isolated, output-only QuanPol gradient instrumentation binary."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


def instrument(source):
    start = source.index("      SUBROUTINE FFOPT2(")
    anchor = "C        - ZERO OFF SOME FORCES -"
    position = source.index(anchor, start)
    position = source.rfind("C\n", start, position)
    dump = """      IF(MASWRK) THEN
         WRITE(IW,'(A,2I10,ES26.17)')
     *   ' QMMM_FORCE_EXPORT_BEGIN ',NAT,NFFAT,ENTOT
         DO IAT=1,NAT
            WRITE(IW,'(A,I10,3ES26.17)') ' QMG ',IAT,
     *      (QMGRD2(III,IAT),III=1,3)
         ENDDO
         DO IFFAT=1,NFFAT
            WRITE(IW,'(A,2I10,3ES26.17)') ' MMG ',IFFAT,
     *      LISTQM(IFFAT),(FFGRD2(III,IFFAT),III=1,3)
         ENDDO
         WRITE(IW,'(A)') ' QMMM_FORCE_EXPORT_END'
      END IF
C
"""
    if "QMMM_FORCE_EXPORT_BEGIN" in source:
        raise ValueError("Source is already instrumented")
    return source[:position] + dump + source[position:]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamess", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    gamess, output = args.gamess.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    objects = output / "object"
    objects.mkdir(exist_ok=True)
    for original in (gamess / "object").iterdir():
        target = objects / original.name
        if original.name != "quanpoc.o" and not target.exists():
            target.symlink_to(original)
    source = gamess / "source/quanpoc.src"
    modified = output / "quanpoc_force.f"
    modified.write_text(instrument(source.read_text()))
    command = ["gfortran", "-c", "-fdefault-integer-8", "-g", "-O2",
               "-fno-aggressive-loop-optimizations", "-funconstrained-commons",
               "-fallow-argument-mismatch", "-I", str(objects), str(modified),
               "-o", str(objects / "quanpoc.o")]
    with (output / "compile.log").open("w") as handle:
        subprocess.run(command, cwd=objects, stdout=handle, stderr=subprocess.STDOUT, check=True)
    config = (gamess / "install.info").read_text()
    config = re.sub(r"(?m)^(setenv GMS_BUILD_DIR\s+).*$", lambda match: match[1] + str(output), config)
    (output / "install.info").write_text(config)
    for name in ("ddi", "tools", "3rd-party", "gms-files.csh", "auxdata", "ddikick.x"):
        if (gamess / name).exists() and not (output / name).exists():
            (output / name).symlink_to(gamess / name, target_is_directory=True)
    with (output / "link.log").open("w") as handle:
        subprocess.run(["tcsh", "-f", str(gamess / "lked"), "gamess", "force"],
                       cwd=output, stdout=handle, stderr=subprocess.STDOUT, check=True)
    launcher = (gamess / "rungms").read_text()
    launcher = launcher.replace(str(gamess), str(output))
    (output / "rungms").write_text(launcher)
    (output / "rungms").chmod(0o755)
    manifest = {"status": "BUILT", "original_source_sha256": sha(source),
                "instrumented_source_sha256": sha(modified),
                "original_binary_sha256": sha(gamess / "gamess.00.x"),
                "binary_sha256": sha(output / "gamess.force.x"),
                "gradient_units": "hartree/bohr", "energy_units": "hartree",
                "instrumentation": "FFOPT2: after QM/MM gradient merge, before masking/fixing",
                "original_installation_modified": False, "compile_command": command}
    (output / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
