:: This script will list the name of all batch files in the gamess-64 folder.
@echo Script                   Description
@echo ======                   ===========
@echo clean-runall-files.bat   Removes exam log, scratch, and restart files.
@echo.
@echo create-parameters.bat    Creates rungms.gms parameter file needed by
@echo                          the rungms.bat batch script.
@echo.
@echo double-click-run.bat     Run all *.inp files from the input folder and
@echo                          store log files in output folder defined in 
@echo                          double-click-run.gms
@echo.
@echo drag-drop-run.bat        Script that accomodates drag-and-drop of *.inp
@echo                          files. Uses [ncpus] and [version] defined in
@echo                          drag-drop-run.gms
@echo.
@echo get-version-names.bat    List available GAMESS binaries and version names
@echo.
@echo runall.bat [version] [ncpus]
@echo                          Run exam files in either serial or parallel
@echo.
@echo rungms.bat [input] [version] [ncpus] [logfile]
@echo                          Run a GAMESS calculation
@echo.
@pause