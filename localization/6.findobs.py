#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:10:10 2024

find the observations near the target pixel

@author: yingyingliu

python 6.find_obs.py 1
"""

import numpy as np
import os

import sys
params_path = '/data42/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case15 as pm
argvs = sys.argv  # command line argument
# Get input&output file from argment
round_ind = argvs[1]

#%% 找出目标网格周围的可用观测数据
if (pm.val_mode()==0):
    filedir = pm.patch_filedir()+'/local_patch/'
    fileout = pm.patch_filedir()+'/obs_patch/'
elif (pm.val_mode()==6):
    filedir = pm.patch_filedir()+'/local_river_patch'+round_ind+'/'
    fileout = pm.patch_filedir()+'/obs_river_patch'+round_ind+'/'
else:
    filedir = pm.patch_filedir()+'/local_select_patch'+round_ind+'/'
    fileout = pm.patch_filedir()+'/obs_select_patch'+round_ind+'/'
os.makedirs(fileout,exist_ok=True)
print(fileout)
indout = pm.patch_filedir()+'/'

def read_txt(fname,i,obsgrid):
    ilon=list()
    ilat=list()
    wt  =list()
    with open(filedir+fname,"rb") as f:
        for line in f:
            variables=list(line.strip().split())
            ilon.append(int(variables[0]))
            ilat.append(int(variables[1]))
            wt.append(float(variables[2]))
    f.close()
    ilon=np.array(ilon,dtype=np.int32)
    ilat=np.array(ilat,dtype=np.int32)
    wt  =np.array(wt,dtype=np.float32)
    obsgrid[i,ilat-1,ilon-1] = wt 
    return obsgrid
 
file_names = os.listdir(filedir)
txt_files = [file for file in file_names if file.endswith('.txt')]

obsgrid = np.zeros((len(txt_files),1320,1500))
countp  = np.zeros((1320,1500),dtype=np.int32)
targtp  = np.zeros((1320,1500),dtype=np.int32)
obs_name = list()
obs_ilon = list()
obs_ilat = list()
i=0
for txt_file in txt_files:
    obs_name.append(txt_file[5:-4])
    obs_ilon.append(txt_file[5:9])
    obs_ilat.append(txt_file[-8:-4])
    obsgrid = read_txt(txt_file,i,obsgrid)
    i=i+1
obs_ilon = np.array(obs_ilon,dtype=np.int32)
obs_ilat = np.array(obs_ilat,dtype=np.int32)

tarlon = list()
tarlat = list()
for i in range(0,1320):
    for j in range(0,1500):
        if np.all(obsgrid[:,i,j]==0.)==True:
            continue
        else:
            countnum = 1
            f_ilat = '{:04d}'.format(i+1)
            f_ilon = '{:04d}'.format(j+1)
            tarlon.append(f_ilon)
            tarlat.append(f_ilat)
            tarname= 'patch'+f_ilon+f_ilat+'.txt' 
            with open(fileout+tarname,'w') as file:         
                for ind in range(0,len(txt_files)):
                    if obsgrid[ind,i,j]==0:
                        continue
                    ind_ilat = '{:04d}'.format(obs_ilat[ind])
                    ind_ilon = '{:04d}'.format(obs_ilon[ind])
                    ind_wt   = '{:.7f}'.format(obsgrid[ind,i,j])
                    countnum = countnum+1
                    if (i==obs_ilat[ind]-1) & (j==obs_ilon[ind]-1) :
                        targtp[i,j]=countnum-1
                    line = f"{ind_ilon}  {ind_ilat}  {ind_wt}\n"
                    file.write(line)
                file.close()
                countnum = countnum-1
            countp[i,j] = countnum

if (pm.val_mode()==0):
    filen=indout+'lon_situ_index.txt'
elif (pm.val_mode()==6):
    filen=indout+'lon_select_river'+round_ind+'.txt'
else:
    filen=indout+'lon_select_index'+round_ind+'.txt'
with open(filen,'w') as file1:
    for row1 in tarlon:
        file1.write(row1 + '\n')

if (pm.val_mode()==0):
    filen=indout+'lat_situ_index.txt'
elif (pm.val_mode()==6):
    filen=indout+'lat_select_river'+round_ind+'.txt'
else:
    filen=indout+'lat_select_index'+round_ind+'.txt'
with open(filen,'w') as file2:
    for row2 in tarlat:
        file2.write(row2 + '\n')

newcount    = np.stack((countp,targtp))
count_flat = newcount.flatten()
with open(fileout+'countnum.bin','wb') as file:
    file.write(count_flat)

