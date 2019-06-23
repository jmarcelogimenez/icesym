call pyrcc5 ../images/images.qrc -o images_rc.py

call dist_win/pyuic ../windows/ICESymMainWindow.ui -o ICESymMainWindow_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "ICESymMainWindow_ui.py"
call dist_win/pyuic ../windows/ConfigurationWidget.ui -o configurationWidget_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "configurationWidget_ui.py"
call dist_win/pyuic ../windows/ValveDialog.ui -o valveDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "valveDialog_ui.py"
call dist_win/pyuic ../windows/AtmosphereDialog.ui -o atmosphereDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "atmosphereDialog_ui.py"
call dist_win/pyuic ../windows/TubeDialog.ui -o tubeDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "tubeDialog_ui.py"
call dist_win/pyuic ../windows/InputDialog.ui -o inputDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "inputDialog_ui.py"
call dist_win/pyuic ../windows/CylinderDialog.ui -o cylinderDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "cylinderDialog_ui.py"
call dist_win/pyuic ../windows/TankDialog.ui -o tankDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "tankDialog_ui.py"
call dist_win/pyuic ../windows/JunctionDialog.ui -o junctionDialog_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "junctionDialog_ui.py"
call dist_win/pyuic ../windows/LogTabWidget.ui -o logTabWidget_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "logTabWidget_ui.py"
call dist_win/pyuic ../windows/PostProcessWidget.ui -o postProcessWidget_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "postProcessWidget_ui.py"
call dist_win/pyuic ../windows/PlotTypeOneWidget.ui -o plotTypeOneWidget_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "plotTypeOneWidget_ui.py"
call dist_win/pyuic ../windows/PlotTypeTwoWidget.ui -o plotTypeTwoWidget_ui.py

%SEDPATH%\sed.exe -i "/setContentsMargins/d" "plotTypeTwoWidget_ui.py"

call dist_win/pyuic ../windows/PlotTypeThreeWidget.ui -o plotTypeThreeWidget_ui.py
%SEDPATH%\sed.exe -i "/setContentsMargins/d" "plotTypeThreeWidget_ui.py"

call dist_win/pyuic ../windows/NewCaseDialog.ui -o newCaseDialog_ui.py
%SEDPATH%\sed.exe -i "/setContentsMargins/d" "newCaseDialog_ui.py"

call dist_win/pyuic ../windows/UsageDialog.ui -o usageDialog_ui.py
%SEDPATH%\sed.exe -i "/setContentsMargins/d" "usageDialog_ui.py"

call dist_win/pyuic ../windows/CurveFormatDialog.ui -o curveFormatDialog_ui.py
%SEDPATH%\sed.exe -i "/setContentsMargins/d" curveFormatDialog_ui.py
