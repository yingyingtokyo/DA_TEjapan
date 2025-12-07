#!/bin/bash

#*** PBS setting when needed
#PBS -q F10
#PBS -l select=1:ncpus=8:mem=160gb
#PBS -j oe
#PBS -m ea
#PBS -V

PBS_O_WORKDIR=/work/a06/yingying/Code/plot/online
cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate data
python 4.plot.py 
