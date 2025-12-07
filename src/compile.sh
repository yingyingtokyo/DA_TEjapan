#!/bin/bash
# compile the online data assimilation code
# ./compile yes to remove the old scripts
##############################
if [ "$1" == "yes" ]; then
    echo "Removing old binaries..."
    rm -rf data_assim make_restart
fi

ifort data_assim.f90 -o data_assim -llapack -lblas
ifort make_restart.f90 -o make_restart
echo "Compilation finished!"
