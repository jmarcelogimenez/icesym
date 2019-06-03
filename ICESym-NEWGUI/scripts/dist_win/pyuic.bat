@Echo Off
If NOT "%ANACONDA3PATH%" == "" (
	@"%ANACONDA3PATH%\python.exe" -m PyQt5.uic.pyuic %1 %2 %3 %4 %5 %6 %7 %8 %9
) else (
	@ECHO The ANACONDA3PATH environment variable was NOT detected.
)