#!/bin/bash

rm -r build/ dist/
pyinstaller --onefile --add-binary dlls/libdl.so.2:. --add-binary dlls/libz.so.1:. --add-binary dlls/libc.so.6:. --add-binary dlls/libGL.so.1:. --add-binary dlls/libxcb-dri3.so.0:. --add-binary dlls/ld-linux-x86-64.so.2:. --icon engine.ico --name icesym ../__main__.py
