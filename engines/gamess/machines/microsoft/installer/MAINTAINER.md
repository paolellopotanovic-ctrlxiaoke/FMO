# Installer Guide

This guide is for the GAMESS maintainer of the Windows 64-bit binary

## Stacksize
Increase stack size for binary
editbin /stack:10000000 C:\cygwin64\home\sarom\gamess\gamess.2023.R1.intel.exe

## Installer
https://www.advancedinstaller.com/

## Simple Package Guide
https://www.advancedinstaller.com/user-guide/tutorial-simple.html

## Notes

New > Installer > Generic : Installer Project : Simple

### Product Details

Name: GAMESS 64-bit 2023 R1 Intel
Version: 2023.1
Publisher: Mark S. Gordon Quantum Theory Group

Support Link: https://github.com/gms-bbg/gamess-issues
Contact: gamess@iastate.edu

### Files and Folders

Application Folder:
- auxdata (folder)
- documentation (folder)
- outputs (folder)
- restart (folder)
- scratch (folder)
- tests (folder)
- windows (folder)
- MS-MPI (folder)
- CITATION.md
- MAINTAINER.md
- RELEASE.md
- gamess.2021.R1.P2.intel.msucc.exe
- double-click-run.bat
- double-click-run.gms
- drag-drop-run.bat
- drag-drop-run.gms
- rungms.bat
- rungms.gms
- Windows-Command-Prompt.lnk
- README.txt
- clean-runall-files.bat
- create-parameters.bat
- get-version-names.bat
- list-contents.bat
- runall.bat

### Install Parameters

Application folder: [PublicFolder]\gamess-64

Package type: 64-bit package for x64 processors

[x] Enable verbose logging

### Builds

[x] Single MSI

MSI name: gamess-64-{name}.msi

### Launch Conditions

[ ] Deselect 32-bit Windows versions

[x] Minimum Physical Memory 1 GB

