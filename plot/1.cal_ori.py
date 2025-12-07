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
import params_ils01 as pm
pm_path = pm.__file__
pm_name = pm.__name__

#%%
def generate_nc(var,type_file):
    dahour = pm.dahour()
    exp_name = pm.runname(pm.mode())
    dahour_str = '{:02d}'.format(dahour) 
    mean_num = pm.ens_mem()
    if 'round' in pm_path:
        savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/round'+str(pm.val_round())+'/'+var+'/'
    elif ('all' in pm_path) & ('case' not in pm_path) :
        savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/all/'+var+'/'
    elif 'case' in pm_path:
        savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/'+pm_name[-5:]+'/'+var+'/all/'
    else:
        savepath = pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/'+var+'/'
    os.makedirs(savepath, exist_ok=True)
    
    start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
    syyyy='%04d' % (start_year)
    smm='%02d' % (start_month)
    sdd='%02d' % (start_date)
    shh='%02d' % (start_hour)
    stime = datetime(start_year,start_month,start_date,start_hour)
    end_year,end_month,end_date,end_hour=pm.endtime() # Start year month date
    eyyyy='%04d' % (end_year)
    emm='%02d' % (end_month)
    edd='%02d' % (end_date)
    ehh='%02d' % (end_hour)
    etime = datetime(end_year,end_month,end_date,end_hour)
    time_step = timedelta(hours=dahour)

    time_len = int((etime-stime).total_seconds() / 3600)
    if 'ils01' in pm_name:
        filepath = pm.DA_dir()+'/'+str(syyyy)+'/out_ils_round'+str(pm.val_round(pm.val_mode()))+'/CaMa_out_'+exp_name+dahour_str+'/'
    else:
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

    def read_bin_off(filepath,ctime,enum):
        ens = '{:03d}'.format(enum)
        fname = filepath+ctime+'_'+ens+'.bin'
        print(fname)
        file_size = os.path.getsize(fname)
        if (file_size>1)&(file_size<8*10**6):
            with open(fname, 'rb') as f:
                data_bin = np.fromfile(f, dtype=np.float32)
                data_re  = np.reshape(data_bin,(dahour,1320,1500))
        else:
            data_re = np.full((dahour,1320,1500),np.nan)
        return data_re

    # assimilated stations
    path_alloc   = pm.alloc_dir()
    if 'round' in pm_path:
        lat_station = '/lat_select_index'+pm_path[-4]+'.txt'
        lat_exclude = '/lat_exclude_index'+pm_path[-4]+'.txt'
        lon_station = '/lon_select_index'+pm_path[-4]+'.txt'
        lon_exclude = '/lon_exclude_index'+pm_path[-4]+'.txt'
        path_off  = pm.outdir() + '/assim_ils' + '{:02d}'.format(pm.dahour()) + '/rivdph/open/rivdph'
    elif 'long' in pm_path:
        lat_station  = '/lat_select_index0.txt'
        lat_exclude  = '/lat_exclude_index0.txt'
        lon_station  = '/lon_select_index0.txt'
        lon_exclude  = '/lon_exclude_index0.txt'
        path_off  = pm.outdir() + '/assim_ils' + '{:02d}'.format(pm.dahour()) + '/rivdph/open/rivdph'
    else:
        lat_station  = '/lat_select_index0.txt'
        lat_exclude  = '/lat_exclude_index0.txt'
        lon_station  = '/lon_select_index0.txt'
        lon_exclude  = '/lon_exclude_index0.txt'
        path_off  = pm.DA_dir()+'/'+str(syyyy)+'/out_ils_round'+str(pm.val_round(pm.val_mode()))+'/assim_ils' + '{:02d}'.format(pm.dahour()) + '/rivdph/open/rivdph' 

    # cross-validated stations
    if 'case' in pm_path:
        lon_point_ind = np.load(pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/lon_case_ind'+pm_name[-2]+pm_name[-1]+'.npy')
        lat_point_ind = np.load(pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/lat_case_ind'+pm_name[-2]+pm_name[-1]+'.npy')
    elif 'round' in pm_path:
        lon_point_ind = np.load(pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/lon_point_round'+pm_path[-4]+'.npy')
        lat_point_ind = np.load(pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/lat_point_round'+pm_path[-4]+'.npy')
    else:
        lon_point_ind = np.load(pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/lon_point_ind.npy')
        lat_point_ind = np.load(pm.plot_dir()+'/exp_'+exp_name+dahour_str+'/lat_point_ind.npy')

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
    data_file = np.full((time_len,len(lat_iy)),np.nan)
    data_off  = np.full((time_len,len(lat_point_ind)),np.nan)
    # only at assimilated stations
    while current_time < etime:
        ctime = current_time.strftime('%Y%m%d%H')
        nyyyy = current_time.strftime('%Y')
        data_day = read_bin(ctime,nyyyy,20,type_file)
        for sta_ind in range(0,len(lat_iy)):
            lat_str = '{:04d}'.format(int(lat_iy[sta_ind]))
            lon_str = '{:04d}'.format(int(lon_ix[sta_ind]))
            data_file[nhour:nhour+dahour,sta_ind] = data_day[:,int(lat_iy[sta_ind])-1,int(lon_ix[sta_ind])-1]
        current_time += time_step
        nhour = nhour+dahour
    # only at validated stations
    nhour = 0
    current_time = stime
    while current_time < etime:
        ctime = current_time.strftime('%Y%m%d%H')
        nyyyy = current_time.strftime('%Y')
        off_day  = read_bin_off(path_off,ctime,20)
        for sta_ind in range(0,len(lat_point_ind)):
            lat_str = '{:04d}'.format(int(lat_point_ind[sta_ind])-1)
            lon_str = '{:04d}'.format(int(lon_point_ind[sta_ind])-1)
            data_off[nhour:nhour+dahour,sta_ind]  = off_day[:,int(lat_point_ind[sta_ind])-1,int(lon_point_ind[sta_ind])-1]
        current_time += time_step
        nhour = nhour+dahour
    data32_ori  = data_file.astype(np.float32)
    off32_ori   = data_off.astype(np.float32)
    mk_bin(data32_ori,'ori')
    mk_bin(off32_ori,'off')

generate_nc('rivdph','C')
#generate_nc('outflw','C')

