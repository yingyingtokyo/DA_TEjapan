#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#*************************************************************************************
# Near real time data assimilation for CaMa-Flood using LETKF and empircal local patches for CaMa-Flood
# [Yingying Liu et al,. (2025)]
# ====================================================================================
# Reference:
# 1. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2020). A framework for estimating 
# global‐scale river discharge by assimilating satellite altimetry. Water Resources Research, 
# 1–34. https://doi.org/10.1029/2020wr027876
# 2. Revel, M., Ikeshima, D., Yamazaki, D., & Kanae, S. (2019). A Physically Based Empirical 
# Localization Method for Assimilating Synthetic SWOT Observations of a Continental-Scale River: 
# A Case Study in the Congo Basin,Water, 11(4), 829. https://doi.org/10.3390/w11040829
# ====================================================================================
# created by Menaka & Yingying Liu
# yingying@iis.u-tokyo.ac.jp
# 2025.10.28
#*************************************************************************************

########################################
#
# This program runs the whole algorithm
#
########################################

## check before run ################
#
# 1. compile ILS
# 2. compile CaMa-Flood
# 3. set parameters params.py@gosh directory
# 4. set local patches folder
# 5. compile fortran codes in ./src folder %sh compile.sh "yes"

import sys
import time
from multiprocessing import Pool
from multiprocessing import Process
import os
#external python codes
import params_case12 as pm
## define of params.py
pm.def_params()
## define of params.py
# ################################################
import wrt_expset as expset
import prep_init as init 
import prep_runoff as inpt
import prep_obs as obs
import main_code as mc

cpu_num = pm.cpu_nums()
os.environ["OMP_NUM_THREADS"] = str(cpu_num)
## write the experimental settings in log file
print ("write experiment log file")
expset.write_text()

## make necessary directories
print ("initializing....")
init.initial()
init.make_initial_infl()

start_time1 = time.time()
## prepare observations
print ("prepare observation")
#obs.download_obs()
obs.prepare_obs()
start_time2 = time.time()
download_time = start_time2 - start_time1
print('download time:'+'{:.4f}'.format(download_time)+'\n')

## initial inflation parameter rho for assimilation
init.make_initial_infl()

## submit ILS scripts to generate ensembles
#############################
# python /data42/yingying/ILS/hagibis/qsub.py
#############################

#print ("preparing ensembles....")
#inpt.ensembles()
#if __name__ == "__main__":  
#    inpt.cal_runoff() 

## prepare localization parameters for DA
# python prep_patch.py 0

## run main code
print ("running main assimilation code....")
mc.main_action()
end_time = time.time()
maincode_time = end_time - start_time2
with open('./experimetal_settings.log','a') as f:
	f.write('download time: '+'{:.4f}'.format(download_time)+'\n')
	f.write('DA time: '+'{:.4f}'.format(maincode_time)+'\n')

## run plot code
print ("HydroDA completed.")
