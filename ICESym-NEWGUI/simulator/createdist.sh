#!/bin/bash

source ../../../bin/activate
rm -r build/ dist/
pyinstaller --onefile --add-binary simCythonCPP.so:. --name exec exec.py
