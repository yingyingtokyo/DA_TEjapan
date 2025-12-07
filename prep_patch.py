#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re
import os
import shutil
import pandas as pd

import sys
params_path = '/data42/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case15 as pm

#*************************************************************************************
# generate empircal local patches for LETKF method
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
# Yingying@IIS 2025
#*************************************************************************************

pm_name = pm.__name__
def def_params():
    pm_line = 'import params_case15 as pm\n'
    def replace_as_pm_line(target_file):
        # replace the definition of pm in target files
        with open(target_file, 'r') as f:
            target_lines = f.readlines()
        new_target_lines = []
        for line in target_lines:
            if 'import params' in line:
                new_target_lines.append(pm_line)
            else:
                new_target_lines.append(line)
        with open(target_file, 'w') as f:
            f.writelines(new_target_lines)
        print("finish definition of params.py！")
    replace_as_pm_line('./localization/4.weight_sfcelv.py')
    replace_as_pm_line('./localization/6.findobs.py')
def_params()

def def_core(num):
    num_str = '{:03d}'.format(num)
    pm_line = "#PBS -l nodes=c"+num_str+":ppn=4 \n"
    def replace_as_pm_line(target_file):
        # replace the definition of pm in target files
        with open(target_file, 'r') as f:
            target_lines = f.readlines()
        new_target_lines = []
        for line in target_lines:
            if 'PBS -l nodes=c' in line:
                new_target_lines.append(pm_line)
            else:
                new_target_lines.append(line)
        with open(target_file, 'w') as f:
            f.writelines(new_target_lines)
        print("finish selecting of core!")
    replace_as_pm_line('./localization/3.make_sft.sh')
    replace_as_pm_line('./localization/4.sub_weight.sh')
    replace_as_pm_line('./localization/5.make_lpara.sh')
    replace_as_pm_line('./localization/6.findobs.sh')
def_core(10)

############ change it according to the needs ############
realmode = pm.real_mode()
realmode = 2 
##########################################################

validatemode = pm.val_mode()
path_alloc   = pm.alloc_dir()
name_alloc   = '/wlv_2019h.xlsx'
df     = pd.read_excel(path_alloc+name_alloc)
print(path_alloc+name_alloc)
df.rename(columns={df.columns[1]: 'latitude', df.columns[2]: 'longitude'}, inplace=True)
id_info =  np.array(df.iloc[:,0])
lat_obs =  np.array(df.iloc[:,1])
lon_obs =  np.array(df.iloc[:,2])
ix     = np.array(df.iloc[:,8])-1  # lon_ind
iy     = np.array(df.iloc[:,9])-1  # lat_ind

argvs = sys.argv
round_ind = argvs[1]
# hagibis: 0-5
# case2: 6-10 

if (validatemode == 2) | (validatemode == 4) | (validatemode==5):
    if (os.path.isfile(path_alloc+'/lat_select_index'+round_ind+'.txt')) & (validatemode!=5):
        print('Dont need regenerate.')
    else:
        fname = pm.patch_filedir()+'/obs_patch/countnum.bin'
        with open(fname, 'rb') as file:
            data_bin = np.fromfile(file, dtype=np.int32)
            data  = np.reshape(data_bin,(2,1320,1500))
            obs_num = data[0,:,:]

        for ind in range(np.shape(ix)[0]):
            if obs_num[iy[ind],ix[ind]] <2:
                lat_obs[ind] = np.nan

        valid_mask = ~np.isnan(lat_obs)
        flood_loc  = pm.flood_reg()
        # typhoon range: 33.4-41.2; 135.4-142
        lat_mask = (lat_obs>flood_loc[0]) & (lat_obs<flood_loc[1])
        lon_mask = (lon_obs>flood_loc[2]) & (lon_obs<flood_loc[3])
        com_mask = valid_mask & lon_mask & lat_mask
        lat_obs_new = lat_obs[com_mask]
        lon_obs_new = lon_obs[com_mask]
        ix_new = ix[com_mask]
        iy_new = iy[com_mask]
        sort_ind = np.argsort(lat_obs_new)
        lat_sort_ind = lat_obs_new[sort_ind]
        lon_sort_ind = lon_obs_new[sort_ind]
        ix_sort_ind  = ix_new[sort_ind]
        iy_sort_ind  = iy_new[sort_ind]

        # reselect stations
        if 'case1' in pm_name:
            num_station = 8
        elif 'case2' in pm_name:
            num_station = 4
        else:
            num_station = 30
        # need to run several times to make sure the number of selected stations is correct
        random_ind = np.zeros(num_station, dtype=int)
        random_num = np.random.choice(len(lat_obs_new), num_station)
        random_ind[:np.shape(random_num)[0]]=random_num
        lat_exclude = iy_sort_ind[random_ind]
        lon_exclude = ix_sort_ind[random_ind]
        print(lat_exclude,lon_exclude)

        mask = ~((np.isin(iy, lat_exclude)) & (np.isin(ix, lon_exclude)))
        lat_filter = iy[mask]
        lon_filter = ix[mask]

        print('DA number: ',np.shape(lat_filter))
        print('validation number: ',np.shape(lat_exclude))

        print(path_alloc+'/lat_select_index'+round_ind+'.txt')
        with open(path_alloc+'/lat_select_index'+round_ind+'.txt', 'w') as file:
            for num in lat_filter:
                file.write(f"{num+1}\n")
            file.close()
        with open(path_alloc+'/lon_select_index'+round_ind+'.txt', 'w') as file:
            for num in lon_filter:
                file.write(f"{num+1}\n")
            file.close()
        # choose for validation
        with open(path_alloc+'/lat_exclude_index'+round_ind+'.txt', 'w') as file:
            for num in lat_exclude:
                file.write(f"{num+1}\n")
            file.close()
        with open(path_alloc+'/lon_exclude_index'+round_ind+'.txt', 'w') as file:
            for num in lon_exclude:
                file.write(f"{num+1}\n")
            file.close()
 
        print("finish selecting!")


patchdir = pm.patch_filedir()
if '06' not in pm_name:
	os.system(f"cp '{path_alloc}/lat_select_index{round_ind}.txt' '{patchdir}/lat_select_index{round_ind}.txt'")
	os.system(f"cp '{path_alloc}/lon_select_index{round_ind}.txt' '{patchdir}/lon_select_index{round_ind}.txt'")
	os.system(f"cp '{path_alloc}/lat_exclude_index{round_ind}.txt' '{patchdir}/lat_exclude_index{round_ind}.txt'")
	os.system(f"cp '{path_alloc}/lon_exclude_index{round_ind}.txt' '{patchdir}/lon_exclude_index{round_ind}.txt'")


if (realmode == 2) | (realmode == 3):
	src = pm.DA_dir()+"/src/"
	patch_size = pm.patch_size()
	outdir = pm.patch_dir()
	thersold = pm.thersold()
	allocdir = pm.alloc_dir()
	camadir  = pm.CaMa_dir()

	os.chdir(src + 'localization/')
	os.system('ifort 1.sfit_sfcelv.f90 -o sfit_sfcelv')

	command1 = (
    f"qsub -v PATCH_SIZE={patch_size},N=696,CAMADIR={camadir},OUTDIR={outdir},ALLOC_DIR={allocdir} "
    "3.make_sft.sh"
).format(
    patch_size=patch_size,
    camadir=camadir,
    outdir=outdir,
    allocdir=allocdir)
	run_job1   = False
	#job_id1   = os.popen(command1).read().strip()
	#run_job1   = True

        # Submit dependent weight job
	if run_job1:
		command2 = f"qsub -W depend=afterok:{job_id1} 4.sub_weight.sh"
	else:
		command2 = f"qsub 4.sub_weight.sh"
	os.system('ifort 5.lpara_sfcelv.f90 -o lpara_sfcelv')
	#job_id2 = os.popen(command2).read().strip()
	run_job2   = False     # True 
	base_qsub = (
    f"qsub "
    f"-v PATCH_SIZE={patch_size},CAMADIR={camadir},OUTDIR={outdir},"
    f"THRESOLD={thersold},ALLOCDIR={allocdir},PATCHDIR={patchdir},"
    f"VALMODE={validatemode},ROUND_IND={round_ind} "
    "5.make_lpara.sh"
)

	if run_job2 and job_id2:
    		command3 = f"qsub -W depend=afterok:{job_id2} " + base_qsub[len("qsub "):]
	else:
		command3 = base_qsub 
	job_id3 = os.popen(command3).read().strip()


	command4 = f"qsub -W depend=afterok:{job_id3} -v ROUND_IND={round_ind} 6.findobs.sh"
	os.system(command4)
