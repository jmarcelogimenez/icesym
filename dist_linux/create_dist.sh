#!/bin/bash

source ../../bin/activate

# check if the enviroment variable is set
if [ -z "$ICESYM_INST_DIR" ]
then
      echo "\$ICESYM_INST_DIR is empty"
      exit
fi

# create the executables
cd $ICESYM_INST_DIR/scripts/dist_linux
./create_dist.sh
cd $ICESYM_INST_DIR/simulator/dist_linux
./create_dist.sh

# create the dist folder
cd $ICESYM_INST_DIR/../dist_linux
rm -r icesym
mkdir icesym
mkdir icesym/simulator
mkdir icesym/runs
mkdir icesym/images
mkdir icesym/etc
cd icesym
cp $ICESYM_INST_DIR/scripts/dist_linux/dist/icesym .
cp $ICESYM_INST_DIR/simulator/dist_linux/dist/exec simulator/.
cp -r $ICESYM_INST_DIR/images/newicons images/.
cp -r $ICESYM_INST_DIR/images/about images/.
cp -r $ICESYM_INST_DIR/etc/bashrc_linux etc/bashrc
cp -r $ICESYM_INST_DIR/cases .
cp -r $ICESYM_INST_DIR/loads .
cp -r $ICESYM_INST_DIR/templates .