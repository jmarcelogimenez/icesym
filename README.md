# ICESym
Authors: Juan Marcelo Gimenez, Ezequiel Jose Lopez, César Pairetti, Santiago Chialvo & Norberto Marcelo Nigro

This project consists of two main parts:

1) The simulator itself, located in the folder ICESym-1D.
2) The GUI, located in the folder ICESym-NEWGUI.

## Linux Distributions

### PACKAGED VERSION

Download this distribution if you are a final user of ICESym, i.e. you will use it only for research or academic propose.

#### Installation

Download the packaged version from [ICESym-Linux64](https://sourceforge.net/projects/icesym/files/latest/download). It includes the simulator and GUI and it was tested in the following distros:

- Ubuntu 14.04
- Ubuntu 16.04
- Ubuntu 18.04

Please, if you use ICESym in another distro please contact us to add it to the list.

#### Usage

For the packaged version, just run:

```bash
source etc/bashrc
./icesym
```

### DEVELOPER VERSION

Follow this steps if you are a developer of ICESym, i.e. you will add new functionalities to the simulator/GUI.

#### Installation

All the following steps are summarized in a bash script that you can download [here](https://sourceforge.net/projects/icesym/files/ICESymInstallScript.sh/download). Note that it seems to be a bug in Ubuntu 18.04 with the default version of PyQt5 used (5.10.1). To avoid this, please use the 5.9.2 version.

Create a virtual environment and download all the needed packages. To start, you will need the following packages:

- Python 2.7
- virtualenv (tested in version 15.1.0)

The first step is to create the python virtual environment, let us call it icesym_virtualenv. The --python flag indicates the path where our python 2.7 is:

```bash
virtualenv --python=/usr/bin/python2.7 icesym_virtualenv
```

The last command created a folder where our virtual enviroment is located. The fortran simulator requires, besides basic python 2.7 packages, numpy, Cython and gfortran compilators. To compile and install the simulator in the folder icesym_virtualenv, follow the next steps:

```bash
export DIR_VIRTUALENV=$HOME/icesym_virtualenv
cd $DIR_VIRTUALENV
source bin/activate
(VENV) git clone https://github.com/jmarcelogimenez/icesym.git
(VENV) pip install numpy==1.8
(VENV) pip install Cython
(VENV) sudo apt-get install gfortran-4.9-multilib
(VENV) cd icesym/ICESym-1D/src
(VENV) make
(VENV) python setup.py install
```

The GUI requires Qt5 (tested with 5.2.1, the default in Ubuntu 14.04), PyQt5 (tested with 5.10.1), sip (tested with 4.19.8), pyinstaller (tested with 3.4) and pyqtgraph. To compile and install all the packages,
follow the next steps (assuming you have already done step 2.2 and have the virtual environment already activated):

*Notes: The NPROC variable depends on the cores in your system that you want to use to compile. The qmake path in the configuration of PyQt5 can vary depending on the distro.
You might need sudo privileges depending where the python libraries are placed.*

```bash
(VENV) sudo apt-get install qt5-default libqt5svg5-dev
(VENV) export NPROCS=4
(VENV) cd $DIR_VIRTUALENV
(VENV) wget -c https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.8/sip-4.19.8.tar.gz
(VENV) tar xvzf sip-4.19.8.tar.gz
(VENV) cd sip-4.19.8/
(VENV) python configure.py --incdir=../include/python2.7
(VENV) make -j $NPROCS
(VENV) make install
(VENV) cd DIR_VIRTUALENV

(VENV) wget -c https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.10.1/PyQt5_gpl-5.10.1.tar.gz
(VENV) tar xvzf PyQt5_gpl-5.10.1.tar.gz
(VENV) cd PyQt-gpl-5.10.1/
(VENV) python configure.py --qmake=/usr/bin/qmake --sip-incdir=../sip-4.19.8/siplib --no-qml-plugin --no-designer-plugin
(VENV) make -j $NPROCS
(VENV) make install
(VENV) cd $DIR_VIRTUALENV

(VENV) pip install -r icesym/dist_linux/requirements.txt
(VENV) cd DIR_VIRTUALENV/icesym/ICESym-NEWGUI/scripts
(VENV) ./compileWidgets_linux.sh
```

There is a final step to use the simulator with the GUI, you must copy the **simCythonCPP.so** library built in ICESym-1D to the simulator/dist_linux folder located in ICESym-NEWGUI. This is done as follows:

```bash
(VENV) cd $DIR_VIRTUALENV/icesym/ICESym-NEWGUI/simulator/dist_linux
(VENV) ./create_dist
```

#### Usage

Running the compiled simulator requires to activate the virtual environment each time:

```bash
cd $DIR_VIRTUALENV
source bin/activate
(VENV) cd icesym/ICESym-NEWGUI
(VENV) python scripts/__main__.py
```

### Creating distributable

To create a linux distributable, assuming you set the virtual enviroment following the steps in Installation (2), just install PyInstaller:

```bash
cd $DIR_VIRTUALENV
source bin/activate
(VENV) pip install PyInstaller==3.4
(VENV) cd $DIR_VIRTUALENV/icesym/ICESym-NEWGUI/
(VENV) source etc/bashrc_linux
(VENV) cd scripts
(VENV) ./compileWidgets_linux.sh
(VENV) cd $DIR_VIRTUALENV/icesym/dist_linux
(VENV) ./create_dist.sh
```

## Windows Distributions

### PACKAGED VERSION

Download this distribution if you are a final user of ICESym, i.e. you will use it only for research or academic propose.

#### Installation

Download the packaged version from [ICESym-Win64](https://sourceforge.net/projects/icesym/files/latest/download) and follow the instructions provided in the text file. Try to locate the root folder within a path that does not contain blank spaces, otherwise the simulator fails to execute inside the GUI. The package was tested in the following distros:

- Windows 7
- Windows 10

Please, if you use ICESym in another distro please contact us to add it to the list.

#### Usage

For packaged version, first execute the etc/bashrc_win.bat and then execute icesym.exe.

### DEVELOPER VERSION

Follow this steps if you are a developer of ICESym, i.e. you will add new functionalities to the simulator/GUI.

#### Installation

Create the work enviroment. First, download and install the needed packages.

- Install the package [Anaconda3-4.3.0-Windows-x86_64](https://repo.continuum.io/archive/). Then, set the ANACONDA3PATH environment variable indicating the installation folder.
- Download and extract [mingw64 4.8](https://sourceforge.net/projects/mingwbuilds/files/host-windows/releases/4.8.0/64-bit/threads-win32/seh/x64-4.8.0-release-win32-seh-rev2.7z/download). Then, add the path of the bin folder inside mingw64 into the PATH environment variable.
- Download and install [Python2.7.10](https://www.python.org/downloads/release/python-2710/). Then, set the PYTHON2PATH environment variable indicating the installation folder.
- Install [sed for Windows](http://gnuwin32.sourceforge.net/packages/sed.htm). Then, set the SEDPATH environment variable of the bin folder path of GnuWin32.

Create the simCythonCPP.pyd file. Inside a cmd console, run:

```bash
%PYTHON2PATH%\python.exe -m pip install numpy
%PYTHON2PATH%\python.exe -m pip install Cython==0.25.2
%PYTHON2PATH%\python.exe setup.py build --compiler=mingw32
```

*Note: If the linking step during the compilation fails, delete the -lmsvcr90 linking library.*

To create the exec.exe file, copy the created .pyd to the folder ICESym-NEWGUI/simulator/dist_win. Then, run:

```bash
%PYTHON2PATH%\python.exe -m pip install PyInstaller==3.4
create_dist.bat
```

#### Usage

For the compiled version:

```bash
cd ICESym-NEWGUI
python scripts/__main__.py
```

### Creating Distributable

To create a windows distributable, assuming you set the work enviroment following the steps in Installation (2), run:

```bash
pip install PyInstaller==3.4
pip install pyqtgraph
conda update setuptools
cd ICESym-NEWGUI/scripts
compileWidgets_windows.bat
```

Then, inside the folder dist_win run:

```bash
create_dist.bat
```
