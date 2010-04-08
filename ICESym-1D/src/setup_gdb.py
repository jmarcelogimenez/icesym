#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext

#!/usr/bin/env python

import os
import numpy
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

#from numpy.distutils.command.build_src import build_src
from numpy.distutils.command import build_src
import Cython
import Cython.Compiler.Main

build_src.Pyrex = Cython
build_src.have_pyrex = True
build_src.build_src.f2py_sources = lambda a,b,c: b

#os.system("cython --cplus --working . -I. simCythonCPP.pyx")

setup(
  name = 'sim_WrapperCythonGDB',
  version      = '0.1.0',
  description  = 'Simulator',
  author       = 'Juan',
  author_email = 'jmarcelogimenez@gmail.com',
  ext_modules=[ 
    Extension(name = "simCythonCPPgdb", 
              sources=[
              #"core.pyx",
              "simCythonCPP.cpp",
              "tube.cc", "junction.cc", "tank.cc", "cylinder.cc", "atmosphere.cc","simulator.cc", "component.cc"
              ,"def_general.f90", "utilities.f90", "gasdyn_utils.f90", "def_simulator.f90", "def_valve.f90", "def_tank.f90",
              "def_junction.f90", "def_cylinder.f90", "def_tube.f90" 
              ], include_dirs = [numpy.get_include()], libraries=['gfortranbegin', 'gfortran'],
			  extra_compile_args=["-g"],extra_link_args=["-g"]
              ),
              
    ],
)
