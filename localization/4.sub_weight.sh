#!/bin/bash
#PBS -l nodes=c010:ppn=4 
#PBS -N local2
#PBS -j oe
#PBS -l walltime=03:30:00
#PBS -V


export PBS_O_WORKDIR=/data50/yingying/HydroDA/src/localization
cd $PBS_O_WORKDIR
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate data
python3 4.weight_sfcelv.py


