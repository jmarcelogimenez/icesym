#!/usr/bin/env bash

set -o errexit # exit when a command fails.
# set -o nounset # exit when your script tries to use undeclared variables.
# set -o xtrace # trace what gets executed. Useful for debugging.

# Set magic variables for current file & dir
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"
__root="$(cd "$(dirname "${__dir}")" && pwd)"
__yellow='\033[0;31m'
__red='\033[0;32m'
__green='\033[0;33m'
__nc='\033[0m' # No Color

usage="$(basename "$0") -- script to instsall ICESym software

where:
    -h  show this help text
    -v  set the folder name of the virtual environment (default: icesym_virtualenv)
    -p  set the python2 PATH (default: /usr/bin/python2.7)
    -n  set the number of processors to compile PyQt5 (default: 1)"

virtualenvFolder=icesym_virtualenv
python2Path=/usr/bin/python2.7
NPROCS=1

while getopts ':hp:v:n:' option; do
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
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

# Creating virtualenv
# ------------------------------------------------------------------------------

sudo apt-get install python2.7-dev
sudo apt-get install python-pip
sudo pip install virtualenv
virtualenv --python=$python2Path $virtualenvFolder

virtualenvPath="${__dir}/$virtualenvFolder"
if [ ! -d "$virtualenvPath" ]; then
  echo -e "${__red}Something wrong happened creating the virtual enviroment.${__nc}"
  exit
else
  echo -e "${__yellow}Virtual enviroment created successfully.${__nc}"
fi

sleep 1

# Compiling fortran simulator
# ------------------------------------------------------------------------------

cd $virtualenvFolder
source "${virtualenvPath}/bin/activate"
git clone https://github.com/jmarcelogimenez/icesym.git
pip install numpy==1.8
pip install Cython
sudo apt-get install gfortran-4.9-multilib
cd icesym/ICESym-1D/src
make
python setup.py install

simCythonPath="${virtualenvPath}/icesym/ICESym-1D/src/build/lib.linux-x86_64-2.7/simCythonCPP.so"
if [ ! -f "$simCythonPath" ]; then
  echo -e "${__red}Something wrong happened compiling the simulator.${__nc}"
  exit
else
  echo -e "${__yellow}Simulator compiled successfully.${__nc}"
fi

sleep 1

# Installing PyQt5
# ------------------------------------------------------------------------------

sudo apt-get install qt5-default libqt5svg5-dev
cd $virtualenvPath
wget -c https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.8/sip-4.19.8.tar.gz
tar xvzf sip-4.19.8.tar.gz
cd sip-4.19.8/
python configure.py --incdir=../include/python2.7
make -j $NPROCS
sudo make install
cd $virtualenvPath

wget -c https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.10.1/PyQt5_gpl-5.10.1.tar.gz
tar xvzf PyQt5_gpl-5.10.1.tar.gz
cd PyQt5_gpl-5.10.1/
python configure.py --qmake=/usr/bin/qmake --sip-incdir=../sip-4.19.8/siplib --no-qml-plugin --no-designer-plugin --confirm-license
make -j $NPROCS
make install
cd $virtualenvPath

echo -e "${__yellow}PyQt5 compiled successfully.${__nc}"

sleep 1

pip install -r icesym/dist_linux/requirements.txt
cd "$virtualenvPath/icesym/ICESym-NEWGUI/scripts"
./compileWidgets_linux.sh

# Creating exec file
# ------------------------------------------------------------------------------

cp $simCythonPath "$virtualenvPath/icesym/ICESym-NEWGUI/simulator/dist_linux"
cd "$virtualenvPath/icesym/ICESym-NEWGUI/simulator/dist_linux"
./create_dist.sh

execPath="${virtualenvPath}/icesym/ICESym-NEWGUI/simulator/dist_linux/dist/exec"
if [ ! -f "$execPath" ]; then
  echo -e "${__red}Something wrong happened creating the exec file.${__nc}"
  exit
else
  echo -e "${__yellow}Exec file created successfully.${__nc}"
fi

cp dist/exec ../.

echo -e "${__green}ICESym installed successfully.${__nc}"
