#!/bin/bash
#*** PBS setting when needed
#PBS -l nodes=c005:ppn=10+c006:ppn=10
#PBS -N qsub_case12
#PBS -l mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -V

echo $LD_LIBRARY_PATH
PBS_O_WORKDIR=/data42/yingying/HydroDA/src
cd $PBS_O_WORKDIR

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate data
cd ./src
ifort data_assim.f90 -o data_assim -L$MKLROOT/lib/intel64_lin -lmkl_rt -mkl
ifort make_restart.f90 -o make_restart

SECONDS=0
/data49/yingying/miniconda/miniconda3/bin/python3.11 run.py
echo "Script execution time: ${SECONDS} seconds" | tee -a time.log
echo "Elapsed time: ${SECONDS} seconds"

