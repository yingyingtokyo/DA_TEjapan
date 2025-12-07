#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:06:49 2024
# calculate the range of ensemble simulation

@author: yingyingliu
"""
#%%
import numpy as np
import netCDF4 as nc
import os
import re
from datetime import datetime, timedelta
from multiprocessing import Pool

import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case14 as pm
pm_path = pm.__file__
pm_name = pm.__name__

#%%
def generate_nc(var,type_file):
    dahour = pm.dahour()
    exp_name = pm.runname(pm.mode())
    dahour_str = '{:02d}'.format(dahour) 
    mean_num = pm.ens_mem()
    if 'all' in pm_path:
        if 'case' in pm_path:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/'+pm_name[-5:]+'/'+var+'/all/'
        else:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/all/'+var+'/'
    else:
        if 'round' in pm_path:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/round'+str(pm.val_round())+'/'+var+'/'
        elif 'case' in pm_path:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/'+pm_name[-6:]+'/'+var+'/valid/'
        elif 'patch' in pm_path:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/patch/'+pm_name[-2:]+'/'+var+'/'
        elif 'ens' in pm_path:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/MEPS/'+var+'/'
        else:
            savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/'+var+'/'
    os.makedirs(savepath, exist_ok=True)
    
    start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
    syyyy='%04d' % (start_year)
    smm='%02d' % (start_month)
    sdd='%02d' % (start_date)
    shh='%02d' % (start_hour)
    if '2022' in syyyy:
        stime = datetime(2022,8,1,0)
    else:
        stime = datetime(start_year,start_month,start_date,start_hour)
    end_year,end_month,end_date,end_hour=pm.endtime() # Start year month date
    eyyyy='%04d' % (end_year)
    emm='%02d' % (end_month)
    edd='%02d' % (end_date)
    ehh='%02d' % (end_hour)
    if '2022' in syyyy:
        etime = datetime(2022,8,6,0)
    else:
        etime = datetime(start_year,start_month,start_date,start_hour)
    time_step = timedelta(hours=dahour)

    time_len = int((etime-stime).total_seconds() / 3600)
    filepath = pm.outdir()+'/CaMa_out_'+exp_name+dahour_str+'/'
    base_date = datetime(1871,1,1,0)
    sta_step = int((stime - base_date).days*24 - 1)
    list_time = [sta_step + t for t in range(time_len)]


    def read_bin(ctime,yyyy,enum,chr_type):
        ens = '{:03d}'.format(enum)
        fname = filepath+ctime+chr_type+ens+'/'+var+yyyy+'.bin'
        print(fname)
        file_size = os.path.getsize(fname)
        if (file_size>1)&(file_size<8*10**6):
            with open(fname, 'rb') as f:
                data_bin = np.fromfile(f, dtype=np.float32)
                data_re  = np.reshape(data_bin,(dahour,1320,1500))
        else:
            data_re = np.full((dahour,1320,1500),np.nan)
        return data_re

    path_alloc   = pm.alloc_dir()
    if 'round' in pm_path:
        lat_station = '/lat_select_index'+pm_path[-4]+'.txt'
        lon_station = '/lon_select_index'+pm_path[-4]+'.txt'
    elif 'long' in pm_path:
        lat_station  = '/lat_select_index0.txt'
        lon_station  = '/lon_select_index0.txt'
    else:
        lat_station  = '/lat_select_index0.txt'
        lon_station  = '/lon_select_index0.txt'
    with open(path_alloc+lat_station, "r", encoding="utf-8") as f:
        lat_iy = f.readlines()
    with open(path_alloc+lon_station, "r", encoding="utf-8") as f:
        lon_ix = f.readlines()

    def mk_bin(data32,typename):
        data_flatten  = data32.flatten()
        file_name = savepath + typename + type_file + '.bin'
        with open(file_name, 'wb') as f:
            f.write(data_flatten.tobytes())
         

    nhour = 0
    current_time = stime
    data_file = np.full((time_len,len(lat_iy),mean_num),np.nan)
    while current_time < etime:
        ctime = current_time.strftime('%Y%m%d%H')
        nyyyy = current_time.strftime('%Y')
        for enum in range(1,1+mean_num):
            data_day = read_bin(ctime,nyyyy,enum,type_file)
            for sta_ind in range(0,len(lat_iy)):
                lat_str = '{:04d}'.format(int(lat_iy[sta_ind]))
                lon_str = '{:04d}'.format(int(lon_ix[sta_ind]))
                data_file[nhour:nhour+dahour,sta_ind,int(enum-1)] = data_day[:,int(lat_iy[sta_ind])-1,int(lon_ix[sta_ind])-1]
        current_time += time_step
        nhour = nhour+dahour
    data32_max  = np.nanmax(data_file,axis=2).astype(np.float32)
    data32_min  = np.nanmin(data_file,axis=2).astype(np.float32)
    data32_mean = np.nanmean(data_file,axis=2).astype(np.float32)
    print(np.shape(data32_mean))
    mk_bin(data32_max,'max')
    mk_bin(data32_min,'min')
    mk_bin(data32_mean,'mean')


generate_nc('rivdph','A')
generate_nc('rivdph','C')
generate_nc('outflw','A')
generate_nc('outflw','C')

