#%%
import numpy as np
import netCDF4 as nc
import os
import re
import pandas as pd
from datetime import datetime, timedelta

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
    filepath = pm.outdir()+'/CaMa_out_'+exp_name+dahour_str+'/'
    base_date = datetime(1871,1,1,0)
    sta_step = int((stime - base_date).days*24 - 1)
    list_time = [sta_step + t for t in range(time_len)]

    def read_bin(ctime,yyyy,enum,chr_type):
        ens = '{:03d}'.format(enum)
        fname = filepath+ctime+chr_type+ens+'/'+var+yyyy+'.bin'
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
        os.makedirs(savepath+'/valid/round'+str(pm.val_round())+'/',exist_ok=True)
        lat_station = '/lat_exclude_index'+pm_path[-4]+'.txt'
        lon_station = '/lon_exclude_index'+pm_path[-4]+'.txt'
    elif 'long' in pm_path:
        os.makedirs(savepath+'/long/',exist_ok=True)
        lat_station  = '/lat_exclude_index0.txt'
        lon_station  = '/lon_exclude_index0.txt'
    elif 'river' in pm_path:
        os.makedirs(savepath+'/river'+str(pm.val_round())+'/',exist_ok=True)
        lat_station  = '/lat_down_index'+str(pm.val_round())+'.txt'
        lat_station1 = '/lat_up_index'+str(pm.val_round())+'.txt'
        lon_station  = '/lon_down_index'+str(pm.val_round())+'.txt'
        lon_station1 = '/lon_up_index'+str(pm.val_round())+'.txt'
    elif 'case' in pm_path:
        lat_station  = '/lat_exclude_index'+str(pm.val_round())+'.txt'
        lon_station  = '/lon_exclude_index'+str(pm.val_round())+'.txt'
    else:
        os.makedirs(savepath+'/valid/',exist_ok=True)
        lat_station  = '/lat_exclude_index0.txt'
        lon_station  = '/lon_exclude_index0.txt'
    with open(path_alloc+lat_station, "r", encoding="utf-8") as f:
        lat_iy = f.readlines()
    with open(path_alloc+lon_station, "r", encoding="utf-8") as f:
        lon_ix = f.readlines()
    if 'river' in pm_path:
        with open(path_alloc+lat_station1, "r", encoding="utf-8") as f:
            lat_iy1 = f.readlines()
        lat_iy = lat_iy1 + lat_iy   # from upstream to downstream
        with open(path_alloc+lon_station1, "r", encoding="utf-8") as f:
            lon_ix1 = f.readlines()
        lon_ix = lon_ix1 + lon_ix   # from upstream to downstream

    data_file = np.full((time_len,len(lat_iy),mean_num),np.nan)

    nhour = 0
    current_time = stime
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

    #os.makedirs(savepath+'/valid/',exist_ok=True)
    for i in range(0,len(lat_iy)):
        data_station = data_file[:,i,:]
        lat_str = '{:04d}'.format(int(lat_iy[i]))
        lon_str = '{:04d}'.format(int(lon_ix[i]))
        for ens in range(0,mean_num):
            ens_str = '{:02d}'.format(ens+1)
            if 'all' in pm_path:
                if 'case' in pm_path:
                    os.makedirs(savepath+'/' + pm_name[-5:] + '/all/',exist_ok=True)
                    file_name = savepath + '/'+ pm_name[-5:]+ '/all/' + lat_str + lon_str + ens_str + type_file +'.bin'            
                else:
                    os.makedirs(savepath+'/all/',exist_ok=True)
                    file_name = savepath + '/all/' + lat_str + lon_str + ens_str + type_file +'.bin'            
            else:
                if 'round' in pm_path:
                    file_name = savepath + '/valid/round'+str(pm.val_round())+'/'+ lat_str + lon_str + ens_str + type_file +'.bin'
                elif 'river' in pm_path:
                    file_name = savepath + '/river'+str(pm.val_round())+'/'+ lat_str + lon_str + ens_str + type_file +'.bin'
                elif 'long' in pm_path:
                    file_name = savepath + '/long/' + lat_str + lon_str + ens_str + type_file +'.bin'
                elif 'case' in pm_path:
                    os.makedirs(savepath+'/'+pm_name[-6:-1]+'/valid'+str(pm.val_round())+'/',exist_ok=True)
                    file_name = savepath + '/'+pm_name[-6:-1]+'/valid'+str(pm.val_round())+'/' + lat_str + lon_str + ens_str + type_file +'.bin'            
                    print(file_name)
                elif 'ens' in pm_path:
                    os.makedirs(savepath+'/ens/',exist_ok=True)
                    file_name = savepath + '/ens/' + lat_str + lon_str + ens_str + type_file +'.bin'            
                else:
                    file_name = savepath + '/valid/'+ lat_str + lon_str + ens_str + type_file +'.bin'
            data32 = data_file[:,i,ens].astype(np.float32)
            data_flatten = data32.flatten()
            with open(file_name, 'wb') as f:
                f.write(data_flatten.tobytes())

generate_nc('rivdph','A')
generate_nc('rivdph','C')
generate_nc('outflw','A')
generate_nc('outflw','C')
