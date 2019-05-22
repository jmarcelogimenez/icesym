#!/bin/bash

rm -r build/ dist/
pyinstaller --onefile --add-binary /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1:. --name icesym __main__.py
