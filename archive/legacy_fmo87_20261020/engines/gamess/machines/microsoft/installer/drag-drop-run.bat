::
::   Mar   08, 2015 :: Sarom Leang :: drag-drop-run.bat
::
:: Allows for the drag and drop of GAMESS input files.
:: Given job.inp the output is written to job.log
::
@IF NOT EXIST drag-drop-run.gms (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: drag-drop-run.gms file not found. Please make one.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Create a new txt file called drag-drop-run.gms in this directory.
  @ECHO This file will contain information about your GAMESS setup.
  @ECHO.
  @ECHO The format of this file is:
  @ECHO VARIABLE=VALUE
  @ECHO.
  @ECHO The file must contain:
  @ECHO   NCPUS=NUMBER OF CPUS
  @ECHO VERSION=GAMESS VERSION
  @ECHO.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Now exiting.
  @EXIT /B
)
@REM
@REM We now set these environmental variables.
@REM
@FOR /F "tokens=1,2 delims==" %%A IN (drag-drop-run.gms) DO @SET %%A=%%B
@REM
@REM Print out what the user entered to invoke this batch file.
@REM

@ECHO -----------------------------------------------------------------------
@ECHO                GAMESS Drag ^& Drop ^& Run Job Submission
@ECHO.
@ECHO                               NCPUS=%NCPUS%
@ECHO                             VERSION=%VERSION%
@ECHO.
@ECHO =======================================================================
@ECHO.

@if [%1]==[] goto :fin
:loop_through_files
@echo rungms.bat %~nx1 %VERSION% %NCPUS% %~n1.log
@call rungms.bat %~nx1 %VERSION% %NCPUS% %~n1.log > NUL
@PING -n 5 127.0.0.1 > NUL
@shift
@if not [%1]==[] goto :loop_through_files
:fin
@echo.
@echo Finish with GAMESS job submission
@echo.
@pause
