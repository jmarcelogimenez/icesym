#!/bin/bash

rm -r build/ dist/
pyinstaller --onefile --add-binary /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1:. --add-binary /usr/lib/x86_64-linux-gnu/libxcb-dri3.so.0:. --icon engine.ico --name icesym ../__main__.py
