#!/bin/bash
# ./sfit_sfcelv 200 720 /work/a06/yingying/CaMa_v411/cmf_v411_pkg/ /work/a06/yingying/out_wlv2019/
#-------------------
#*** PBS setting when needed
#PBS -l nodes=c010:ppn=4 
#PBS -N local1
#PBS -j oe
#PBS -m ea
#PBS -V
export PBS_O_WORKDIR=/data50/yingying/HydroDA/src/localization
cd $PBS_O_WORKDIR
#-----------------

./sfit_sfcelv $PATCH_SIZE $N $CAMADIR $OUTDIR $ALLOC_DIR 




