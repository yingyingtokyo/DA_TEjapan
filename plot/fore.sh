#!/bin/bash

#*** PBS setting when needed
#PBS -l nodes=c002:ppn=10+c003:ppn=10
#PBS -l mem=160gb
#PBS -j oe
#PBS -m ea
#PBS -V

source ~/.bashrc
PBS_O_WORKDIR=/data42/yingying/HydroDA/src/plot
cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate data
python 7.make_24case.py 
