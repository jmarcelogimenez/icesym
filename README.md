# ICESym
Authors: Juan Marcelo Gimenez, Ezequiel Jose Lopez, Santiago Chialvo & Norberto Marcelo Nigro

This project consist in two main parts:

1) The simulator itself, located in the folder ICESym-1D.
2) The GUI, located in the folder ICESym-NEWGUI.

## Installation

There are two ways to install ICESym:

1) Download the packaged version from: [ICESym](https://sourceforge.net/projects/icesym/). It includes the simulator and GUI and it was tested in the following distros:

- Ubuntu 14.04 
- Ubuntu 16.04
- Ubuntu 18.04

Please, if you use ICESym in another distro please contact us to add it to the list.

2) Create a virtual enviroment and download all the needed packages. To start, you will need the following packages:

- Python 2.7
- virtualenv (tested in version 15.1.0)

2.1) The first step is to create the python virtual enviroment, let us call it DIR_VIRTUALENV. The --python flag indicates the path where our python 2.7 is:

```bash
virtualenv --python=/usr/bin/python2.7 DIR_VIRTUALENV
```

2.2) The simulator requires, besides basic python 2.7 packages, numpy, Cython and gfortran compilators. To compile and install, follow the next steps:

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

2.3) The GUI requires Qt5 (tested with 5.2.1, default in Ubuntu 14.04), PyQt5 (tested with 5.10.1), sip (tested with 4.19.8), pyinstaller (tested with 3.4) and pyqtgraph. To compile and install all the packages, 
follow the next steps (assuming you have already done step 2.2 and have the virtual enviroment already activated):

*Notes: The NPROC variable depends on the cores in your system that you want to use to compile. The qmake path in the configuration of PyQt5 can vary depending on the distro.
You might need sudo privileges depending where the python libraries are placed.*

```bash
sudo apt-get install qt5-default
cd DIR_VIRTUALENV
wget -c https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.8/sip-4.19.8.tar.gz
tar xvzf sip-4.19.8.tar.gz 
cd sip-4.19.8/
python configure.py --incdir=../../include/python2.7
make -jNPROC
make install

wget -c https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.10.1/PyQt5_gpl-5.10.1.tar.gz
tar xvzf PyQt-gpl-5.10.1.tar.gz 
cd PyQt-gpl-5.10.1/
python configure.py --qmake=/usr/bin/qmake --sip-incdir=../sip-4.19.8/siplib --no-qml-plugin --no-designer-plugin
make -jNPROC
make install

pip install pyqtgraph
pip install pyinstaller==3.4
```

There is a final step to use the simulator with the GUI, you must copy the **simCythonCPP.so** library builded in ICESym-1D to the simulator folder located in ICESym-NEWGUI. Then:

```bash
cd ICESym-NEWGUI/simulator
./createdist
```

## Usage

1) For packaged version, just run:

```bash
./icesym
```

2) For the compiled version:

```bash
cd ICESym-NEWGUI
python scripts/__main__.py
```