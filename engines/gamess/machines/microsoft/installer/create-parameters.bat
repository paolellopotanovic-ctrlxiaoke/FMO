::
::   Mar   08, 2015 :: Sarom Leang :: created create-parameters.bat
::
:: Create rungms.gms parameter file
::
@for %%f in ("%cd%") do @(
@echo GAMESSDIR=%%~f>rungms.gms
@echo AUXDATADIR=%%~f\auxdata>>rungms.gms
@echo RESTARTDIR=%%~f\restart>>rungms.gms
@echo SCRATCHDIR=%%~f\scratch>>rungms.gms
@echo GAMESSDIR=%%~f
@echo AUXDATADIR=%%~f\auxdata
@echo RESTARTDIR=%%~f\restart
@echo SCRATCHDIR=%%~f\scratch)
@echo.
@echo Finished rungms.gms parameter file
@echo.
@pause