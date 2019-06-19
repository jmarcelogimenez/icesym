# ICESym
Authors: Juan Marcelo Gimenez, Ezequiel Jose Lopez, César Pairetti, Santiago Chialvo & Norberto Marcelo Nigro

This project consists of two main parts:

1) The simulator itself, located in the folder ICESym-1D.
2) The GUI, located in the folder ICESym-NEWGUI.

## Linux Distributions

### PACKAGED VERSION:

#### Installation

Download the packaged version from [ICESym-Linux64](https://sourceforge.net/projects/icesym/files/icesym-l64-20190524.tar.gz/download). It includes the simulator and GUI and it was tested in the following distros:

- Ubuntu 14.04 
- Ubuntu 16.04
- Ubuntu 18.04

Please, if you use ICESym in another distro please contact us to add it to the list.

#### Usage

For packaged version, just run:

```bash
source etc/bashrc
./icesym
```

### FROM SOURCES VERSION:

#### Installation

Create a virtual environment and download all the needed packages. To start, you will need the following packages:

- Python 2.7
- virtualenv (tested in version 15.1.0)

The first step is to create the python virtual environment, let us call it DIR_VIRTUALENV. The --python flag indicates the path where our python 2.7 is:

```bash
virtualenv --python=/usr/bin/python2.7 DIR_VIRTUALENV
```

The simulator requires, besides basic python 2.7 packages, numpy, Cython and gfortran compilators. To compile and install, follow the next steps:

```bash
cd DIR_VIRTUALENV
source bin/activate
git clone https://github.com/jmarcelogimenez/icesym.git
pip install numpy==1.8
pip install Cython
sudo apt-get install gfortran-4.9-multilib
cd icesym/ICESym-1D/src
make
python setup.py install
```

The GUI requires Qt5 (tested with 5.2.1, the default in Ubuntu 14.04), PyQt5 (tested with 5.10.1), sip (tested with 4.19.8), pyinstaller (tested with 3.4) and pyqtgraph. To compile and install all the packages, 
follow the next steps (assuming you have already done step 2.2 and have the virtual environment already activated):

*Notes: The NPROC variable depends on the cores in your system that you want to use to compile. The qmake path in the configuration of PyQt5 can vary depending on the distro.
You might need sudo privileges depending where the python libraries are placed.*

```bash
sudo apt-get install qt5-default
cd DIR_VIRTUALENV
wget -c https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.8/sip-4.19.8.tar.gz
tar xvzf sip-4.19.8.tar.gz 
cd sip-4.19.8/
python configure.py --incdir=../include/python2.7
make -j
make install
cd DIR_VIRTUALENV

wget -c https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.10.1/PyQt5_gpl-5.10.1.tar.gz
tar xvzf PyQt5_gpl-5.10.1.tar.gz
cd PyQt-gpl-5.10.1/
python configure.py --qmake=/usr/bin/qmake --sip-incdir=../sip-4.19.8/siplib --no-qml-plugin --no-designer-plugin
make -j
make install
cd DIR_VIRTUALENV

pip install -r icesym/dist_linux/requirements.txt
cd DIR_VIRTUALENV/icesym/ICESym-NEWGUI/scripts
./compileWidgets_linux.sh
```

There is a final step to use the simulator with the GUI, you must copy the **simCythonCPP.so** library built in ICESym-1D to the simulator/dist_linux folder located in ICESym-NEWGUI. Then:

```bash
cd ICESym-NEWGUI/simulator/dist_linux
./create_dist
```

#### Usage

For the compiled version:

```bash
cd ICESym-NEWGUI
python scripts/__main__.py
```

### Creating distributable

To create a linux distributable, assuming you set the virtual enviroment following the steps in Installation (2), just install PyInstaller:

```bash
pip install PyInstaller==3.4
```

Then, inside the folder dist_linux:

```bash
./create_dist.sh
```

## Windows Distributions

### PACKAGED VERSION

#### Installation

Download the packaged version from [ICESym-Win64](https://sourceforge.net/projects/icesym/files/icesym-w64-20190528.zip/download). It includes the simulator and GUI and it was tested in the following distros:

- Windows 7
- Windows 10

Please, if you use ICESym in another distro please contact us to add it to the list.

#### Usage

For packaged version, first execute the etc/bashrc_win.bat and then execute icesym.exe.

### FROM SOURCES VERSION:

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
