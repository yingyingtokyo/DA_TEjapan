#!/bin/bash

#*** PBS setting when needed
#PBS -l nodes=c003:ppn=2+c005:ppn=2
#PBS -l mem=160gb
#PBS -j oe
#PBS -m ea
#PBS -V

source ~/.bashrc
PBS_O_WORKDIR=/data42/yingying/HydroDA/src/plot
cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate data
python 1.cal_ens_range.py 
