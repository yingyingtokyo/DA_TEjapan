#!/bin/bash

#*** PBS setting when needed
#PBS -l nodes=node37:ppn=8
#PBS -j oe
#PBS -m ea
#PBS -V

PBS_O_WORKDIR=/work/a06/yingying/Code/plot/online
cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate data
python 4.plot_offline.py 
