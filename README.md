# ICESym
Automatically exported from code.google.com/p/icesym
Authors: Juan Marcelo Gimenez, Ezequiel Jose Lopez & Norberto Marcelo Nigro

## Compilation

ICESym requires, besides basic python 2.7 packages, gcc and gfortran compilators. To compile, follow the next steps:

~~bash~~
cd ICESym-1D/src
make
python setup.py install
~~~~

You might need sudo privileges depending where the python libraries are placed.

## Documentation
The template,py file in *ICESym-1D/src/* folder has a brief description of each parameter for each element.

## Tutorials
A simulation can be generated by Graphical User Interface (GUI) or manually with a text editor. In both cases, the setup is defined as a python dictionary in a .py file.

### Graphical User Interface tutorials
The GUI can be started by running:
~~bash~~
cd ICESym-GUI
python Simulator GUI.py
~~~~

The following GUI generated tutorials are placed at *ICESym-GUI/saves/*:

*QUB400.py:
*multicyl4TSI.py:
*monocyl4TSI.py:
*compBiCyl.py:
*KamAZ7405.py:

To **run** any of these cases, open it with *File -> Open* option and then go to *Simulation -> Run* option.

### Manually generated tutorials
For parametrization and optimization purposes, it is sometimes more convenient to edit the python file directly.
Text editor tutorials are placed at *ICESym-1D/src/*:

*compBicCyl.py:
*test6.py:
*opse.py:
*ctest6.py:
*BoscaMono.py:

To run any of these, just follow the next lines:
~~bash~~
cd ICESym-1D/src
python main.py
~~~~

The prompt will then ask for the case name.