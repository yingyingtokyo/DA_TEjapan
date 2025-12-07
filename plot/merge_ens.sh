#!/bin/bash

#*** PBS setting when needed
#PBS -q E10
#PBS -l select=1:ncpus=4:mem=60gb
#PBS -j oe
#PBS -m ea
#PBS -V

PBS_O_WORKDIR=/work/a06/yingying/Code/plot/online
cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate data
python 2.merge_ens.py 
