#!/bin/bash
#-------------------
#*** PBS setting when needed
#PBS -l nodes=c010:ppn=4 
#PBS -N local3
#PBS -j oe
#PBS -l walltime=01:30:00
#PBS -V
echo $LD_LIBRARY_PATH
export PBS_O_WORKDIR=/data50/yingying/HydroDA/src/localization
cd $PBS_O_WORKDIR

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate data

echo "valmode:"$VALMODE
gfortran 5.lpara_sfcelv.f90 -o lpara_sfcelv

if (( $(echo "$THRESOLD > 0.6" | bc -l) )) && (( $(echo "$THRESOLD <= 0.7" | bc -l) )); then
    GAUSS_THRESOLD=0.1
elif (( $(echo "$THRESOLD > 0.7" | bc -l) )) && (( $(echo "$THRESOLD <= 0.8" | bc -l) )); then
    GAUSS_THRESOLD=0.25
elif (( $(echo "$THRESOLD > 0.8" | bc -l) )) && (( $(echo "$THRESOLD <= 0.9" | bc -l) )); then
    GAUSS_THRESOLD=0.4
else
    GAUSS_THRESOLD=0.0
fi

echo $PATCH_SIZE
echo $CAMADIR
echo $OUTDIR
echo $THRESOLD
echo $GAUSS_THRESOLD
echo $ALLOCDIR
echo $PATCHDIR
echo $VALMODE
echo $ROUND_IND

./lpara_sfcelv "$PATCH_SIZE" "$CAMADIR" "$OUTDIR" "$THRESOLD" "$GAUSS_THRESOLD" "$ALLOCDIR" "$PATCHDIR" "$VALMODE" "$ROUND_IND" 



