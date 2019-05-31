@Echo Off
If "%PYTHON2PATH%" == "" GOTO NOPATH
:YESPATH
	@ECHO The PYTHON2PATH environment variable was detected.
	%PYTHON2PATH%\python.exe -m PyInstaller --onefile --name exec exec.spec
GOTO :END
:NOPATH
	@ECHO The PYTHON2PATH environment variable was NOT detected.
GOTO :END
:END