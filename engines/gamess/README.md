# General ab initio Quantum Chemistry Software Package (GAMESS)
# http://www.msg.chem.iastate.edu/

## Quick Start (Public Release)

Quick start to building GAMESS:

1.  Generate an **install.info** file for your system. Two approaches are available:
   - (Method 1) Run `./config` and provide answers to the questions asked about your system configuration. Always run config, copy-pasting an old install.info doesn't work! 
   - (Method 2) Run `./bin/create-install-info.py` and provide appropriate command-line flags. 
2.  `make ddi`
3.  (Optional) if LibXC interface is enabled
    - `make libxc -j$(nproc)`
4.  (Optional) if your math library does not provide its own LAPACK routines (e.g., not LAPACK complete)
    - `./tools/lapack/download-lapack.csh`
    - `make lapack`
5.  `make`
    - You can invoke `make` in parallel to speed up the build process by invoking `make -j`

The config script will ask you for the paths for the needed tools for compilation, this include:

 - Math library
 - MPI install

## GAMESS software documentation
* Section 1 - docs-intro.txt - Overview
* Section 2 - docs-input.txt - Input Description
* Section 3 - docs-tests.txt - Input Examples
* Section 4 - docs-references.txt  - Further Information
* Section 5 - docs-prog.txt  - Programmer's Reference
* Section 6 - docs-hardware.txt  - Hardware Specifics

## Quick Start (Development Branch)

How to:
* [get the latest development source code](https://github.com/gms-bbg/gamess/wiki/Getting-the-Latest-Development-Source-Code).
* [checkout a branch](https://github.com/gms-bbg/gamess/wiki/Checking-out-a-branch).
* [submit a new feature to GAMESS](https://github.com/gms-bbg/gamess/wiki/Submitting-a-New-Feature-to-GAMESS).
* [submit a bug-report](https://github.com/gms-bbg/gamess/wiki/Submitting-a-Bug-Report).

Developer resources:
* [wiki](https://github.com/gms-bbg/gamess/wiki)
* [coding policy](https://github.com/gms-bbg/gamess/blob/development/DEVELOPERS.md)

## How to cite GAMESS

* See [CITATION.md](CITATION.md)
