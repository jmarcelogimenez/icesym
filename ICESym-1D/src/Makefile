.PHONY: default src build install clean

default: build

src: simCythonCPP.cpp

simCythonCPP.cpp: simCythonCPP.pyx $(wildcard *.pyx)
	cython --cplus --working . -I. $<

build: src
	python setup.py build

install: build
	python setup.py install

install_home: build
	python setup.py install --home=${HOME}

clean:
	python setup.py clean --all
	${RM} simCythonCPP.cpp *.pyc *~
