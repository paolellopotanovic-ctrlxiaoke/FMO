::
::   Mar   08, 2015 :: Sarom Leang :: get-version-names.bat
::
:: Looks in current directory for gamess.*.exe and outputs the
:: version names to be used in 'rungms.bat' or 'rungms-serial.bat'
::
@setlocal enabledelayedexpansion
@echo ------					-------
@echo Binary					Version
@echo ======					=======
@for %%f in (gamess.*.exe) do @(
@set temp=%%~nf
@set version=!temp:gamess.=!
@echo %%f		!version!)
@echo.
@pause
