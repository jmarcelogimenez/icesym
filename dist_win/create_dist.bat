rd /s /q icesym
mkdir icesym
mkdir icesym\cases
xcopy ..\ICESym-NEWGUI\cases icesym\cases
mkdir icesym\images
mkdir icesym\images\newicons
mkdir icesym\images\about
xcopy ..\ICESym-NEWGUI\images\newicons icesym\images\newicons
xcopy ..\ICESym-NEWGUI\images\about icesym\images\about
mkdir icesym\loads
xcopy ..\ICESym-NEWGUI\loads icesym\loads
mkdir icesym\docs
xcopy ..\ICESym-NEWGUI\docs icesym\docs
mkdir icesym\runs
mkdir icesym\simulator
rd /s /q ..\ICESym-NEWGUI\simulator\dist_win\dist
rd /s /q ..\ICESym-NEWGUI\simulator\dist_win\build
cd ..\ICESym-NEWGUI\simulator\dist_win\
call create_dist.bat
cd %~dp0%
xcopy ..\ICESym-NEWGUI\simulator\dist_win\dist\exec.exe icesym\simulator
mkdir icesym\templates
xcopy ..\ICESym-NEWGUI\templates icesym\templates
rd /s /q ..\ICESym-NEWGUI\scripts\dist_win\dist
rd /s /q ..\ICESym-NEWGUI\scripts\dist_win\build
cd ..\ICESym-NEWGUI\scripts\dist_win\
call create_dist.bat
cd %~dp0%
xcopy ..\ICESym-NEWGUI\scripts\dist_win\dist\icesym.exe icesym
xcopy ..\ICESym-NEWGUI\etc\bashrc_win.bat icesym
pause
