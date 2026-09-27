::
::   Mar   08, 2015 :: Sarom Leang :: created clean-runall-files.bat
::
:: Cleans us log files, restart files, and scratch files produced
:: from a 'runall.bat' or 'runall-serial.bat' 
::
@erase /Q exam*.log 2>NUL
@erase /Q tests\travis-ci\exam*.log 2>NUL
@erase /Q tests\travis-ci\parallel\exam*.log 2>NUL
@erase /Q scratch\exam*.* 2>NUL
@erase /Q restart\exam*.* 2>NUL
@echo.
@echo Finish removing exam*log, scratch\exam*.*, and restart\exam*.* files
@echo.
@pause
