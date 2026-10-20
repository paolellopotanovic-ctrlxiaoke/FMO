::
::   Mar   08, 2015 :: Sarom Leang :: double-click-run.bat
::
:: Run all *.inp files in directory INPUT_DIRECTORY and writes them to
:: *.log files in directory OUTPUT_DIRECTORY
::
@IF NOT EXIST double-click-run.gms (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: double-click-run.gms file not found. Please make one.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Create a new txt file called drag-drop-run.gms in this directory.
  @ECHO This file will contain information about your GAMESS setup.
  @ECHO.
  @ECHO The format of this file is:
  @ECHO VARIABLE=VALUE
  @ECHO.
  @ECHO The file must contain:
  @ECHO            NCPUS=NUMBER OF CPUS
  @ECHO          VERSION=GAMESS VERSION
  @ECHO  INPUT_DIRECTORY=LOCAL DIRECTORY NAME CONTAINING INPUT FILES
  @ECHO OUTPUT_DIRECTORY=LOCAL DIRECTORY NAME TO PUT OUTPUT FILES
  @ECHO.
  @ECHO INPUT_DIRECTORY AND OUTPUT_DIRECTORY MUST EXIST WITHIN THE CURRENT FOLDER
  @ECHO.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Now exiting.
  @EXIT /B
)
@REM
@REM We now set these environmental variables.
@REM
@FOR /F "tokens=1,2 delims==" %%A IN (double-click-run.gms) DO @SET %%A=%%B
@REM
@REM Print out what the user entered to invoke this batch file.
@REM

@ECHO -----------------------------------------------------------------------
@ECHO                 GAMESS Double Click ^& Run Job Submission
@ECHO.
@ECHO                               NCPUS=%NCPUS%
@ECHO                             VERSION=%VERSION%
@ECHO                     INPUT_DIRECTORY=%INPUT_DIRECTORY%
@ECHO                    OUTPUT_DIRECTORY=%OUTPUT_DIRECTORY%
@ECHO.
@ECHO =======================================================================
@ECHO.

@for /f %%f in ('dir /b %INPUT_DIRECTORY%\*.inp') do @(
@copy %INPUT_DIRECTORY%\%%f . >NUL
@echo rungms.bat %%f %VERSION% %NCPUS% ^> %OUTPUT_DIRECTORY%\%%~nf.log 2^>NUL
@call rungms.bat %%f %VERSION% %NCPUS% > %OUTPUT_DIRECTORY%\%%~nf.log 2>NUL
@PING -n 5 127.0.0.1 >NUL 2>NUL
@erase /Q %%f >NUL 2>NUL)
@goto fin

:fin
@echo.
@echo Finish with GAMESS job submission
@echo.
@pause
