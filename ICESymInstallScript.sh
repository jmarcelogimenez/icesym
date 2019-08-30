#!/usr/bin/env bash


# Useful functions and variables
# ------------------------------------------------------------------------------

# Check for specific package
# Return values:
#  0 - package is installed
#  1 - package is not installed, it is available in package repository
#  2 - package is not installed, it is not available in package repository
check_for_package() {
  if dpkg-query -s "$1" 1>/dev/null 2>&1; then
    return 0   # package is installed
  else
    if apt-cache show "$1" 1>/dev/null 2>&1; then
      return 1 # package is not installed, it is available in package repository
    #else
      #return 2 # package is not installed, it is not available in package repository
    fi
  fi
}

set -o errexit # exit when a command fails.
# set -o nounset # exit when your script tries to use undeclared variables.
# set -o xtrace # trace what gets executed. Useful for debugging.

# Set magic variables for current file & dir
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"
__root="$(cd "$(dirname "${__dir}")" && pwd)"
__red='\033[0;31m'
__green='\033[0;32m'
__yellow='\033[0;33m'
__nc='\033[0m' # No Color

# Script usage information
# ------------------------------------------------------------------------------

usage="$(basename "$0") -- script to instsall ICESym software

where:
    -h  show this help text
    -v  set the folder name of the virtual environment (default: icesym_virtualenv)
    -p  set the python2 PATH (default: /usr/bin/python2.7)
    -q  set the version of PyQt5 (default: 5.10.1)
    -n  set the number of processors to compile PyQt5 (default: 1)"

virtualenvFolder=icesym_virtualenv
python2Path=/usr/bin/python2.7
NPROCS=1
PYQT5_VERSION=5.10.1

while getopts ':hp:v:n:q:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    p) python2Path=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    v) virtualenvFolder=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    n) NPROCS=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    q) PYQT5_VERSION=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

# Installing all needed packages
# ------------------------------------------------------------------------------

packages="python2.7-dev python-pip gfortran-5-multilib qt5-default libqt5svg5-dev"

for package in $packages; do
  if check_for_package "$package"; then
    echo -e "${__green} "$package" package is installed, skipping.${__nc}"
  else
    echo -e "${__yellow} "$package" package is not installed, installing.${__nc}"
    sudo apt-get install $package
  fi
done

# Creating virtualenv
# ------------------------------------------------------------------------------

sudo pip install virtualenv
virtualenvPath="${__dir}/$virtualenvFolder"
createVirtualEnv=true
if [ -d "$virtualenvPath" ]; then
  read -p "There is already a folder in that location. Do you want to erase it and install all again? (y/n) " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    sudo rm -r $virtualenvPath
  else
    createVirtualEnv=false
  fi
fi

if [ "$createVirtualEnv" = true ]; then
  virtualenv --python=$python2Path $virtualenvFolder

  if [ ! -d "$virtualenvPath" ]; then
    echo -e "${__red}Something wrong happened creating the virtual enviroment.${__nc}"
    exit
  else
    echo -e "${__green}Virtual enviroment created successfully.${__nc}"
  fi
else
  echo -e "${__yellow}Virtual enviroment already created, skipping.${__nc}"
fi

sleep 1

# Compiling fortran simulator
# ------------------------------------------------------------------------------

cd $virtualenvFolder
source "${virtualenvPath}/bin/activate"

ICESymPath="${__dir}/$virtualenvFolder/icesym"
cloneICESym=true
if [ -d "$ICESymPath" ]; then
  read -p "There is already a repository cloned inside the virtual environment folder. Do you want to erase it and install all again? (y/n) " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    sudo rm -r $ICESymPath
  else
    cloneICESym=false
  fi
fi

if [ "$cloneICESym" = true ]; then
  git clone https://github.com/jmarcelogimenez/icesym.git
fi

pip install numpy==1.8
pip install Cython

simCythonPath="${virtualenvPath}/icesym/ICESym-1D/src/build/lib.linux-x86_64-2.7/simCythonCPP.so"
compileSimulator=true
if [[ -f "$simCythonPath" ]]; then
  read -p "There is already a compiled version of the simulator. Do you want to erase it and compile it again? (y/n) " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    compileSimulator=false
  fi
fi

if [ "$compileSimulator" = true ]; then
  gfortanSymLink=/usr/bin/gfortran
  libgfortranSymLink=/usr/lib/libgfortran.so

  if [[ ! -f "$gfortanSymLink" ]]; then
    sudo ln -s /usr/bin/gfortran-5 /usr/bin/gfortran
  fi
  if [[ ! -f "$libgfortranSymLink" ]]; then
    sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/5/libgfortran.so /usr/lib/libgfortran.so
  fi

  cd icesym/ICESym-1D/src
  make
  python setup.py install

  if [ ! -f "$simCythonPath" ]; then
    echo -e "${__red}Something wrong happened compiling the simulator.${__nc}"
    exit
  else
    echo -e "${__green}Simulator compiled successfully.${__nc}"
  fi
fi
sleep 1

# Installing PyQt5
# ------------------------------------------------------------------------------

cd $virtualenvPath

sipPath="${virtualenvPath}/sip-4.19.8.tar.gz"
if [[ ! -f "$sipPath" ]]; then
  wget -c https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.8/sip-4.19.8.tar.gz
  tar xvzf sip-4.19.8.tar.gz
  cd sip-4.19.8/
  python configure.py --incdir=../include/python2.7
  make -j $NPROCS
  sudo make install
else
  echo -e "${__yellow}A sip tar has been found. If you want to re-compile sip, please erase or move it.${__nc}"
fi

cd $virtualenvPath
# 5.9.2
PyQt5Path="${virtualenvPath}/PyQt5_gpl-"$PYQT5_VERSION".tar.gz"
if [[ ! -f "$PyQt5Path" ]]; then
  wget -c https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-$PYQT5_VERSION/PyQt5_gpl-$PYQT5_VERSION.tar.gz
  tar xvzf PyQt5_gpl-$PYQT5_VERSION.tar.gz
  cd PyQt5_gpl-$PYQT5_VERSION/
  python configure.py --qmake=/usr/bin/qmake --sip-incdir=../sip-4.19.8/siplib --no-qml-plugin --no-designer-plugin --confirm-license
  make -j $NPROCS
  make install
  echo -e "${__green}PyQt5 compiled successfully.${__nc}"
else
  echo -e "${__yellow}A PyQt5 tar has been found. If you want to re-compile PyQt5, please erase or move it.${__nc}"
fi

cd $virtualenvPath
sleep 1

pip install -r icesym/dist_linux/requirements.txt
cd "$virtualenvPath/icesym/ICESym-NEWGUI/scripts"
./compileWidgets_linux.sh

# Creating exec file
# ------------------------------------------------------------------------------

execPath="${virtualenvPath}/icesym/ICESym-NEWGUI/simulator/dist_linux/dist/exec"
createExec=true
if [[ -f "$execPath" ]]; then
  read -p "There is already an exec created. Do you want to erase it and create it again? (y/n) " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    createExec=false
  fi
fi

if [ "$createExec" = true ]; then
  cp $simCythonPath "$virtualenvPath/icesym/ICESym-NEWGUI/simulator/dist_linux"
  cd "$virtualenvPath/icesym/ICESym-NEWGUI/simulator/dist_linux"
  ./create_dist.sh

  execPath="${virtualenvPath}/icesym/ICESym-NEWGUI/simulator/dist_linux/dist/exec"
  if [ ! -f "$execPath" ]; then
    echo -e "${__red}Something wrong happened creating the exec file.${__nc}"
    exit
  else
    echo -e "${__green}Exec file created successfully.${__nc}"
  fi
  cp dist/exec ../.
fi

echo -e "${__green}ICESym installed successfully.${__nc}"
