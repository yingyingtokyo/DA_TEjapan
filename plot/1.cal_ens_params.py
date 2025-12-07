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

import sys
params_path = '/data42/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params as pm

#%%
def generate_nc(var,type_file):
    dahour_str = '{:02d}'.format(dahour) 
    mean_num = 20
    savepath = '/work/a06/yingying/Code/plot/exp_'+exp_name+dahour_str+'/'+var+'/'
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
    data_file = np.full((time_len,1320,1500,mean_num),np.nan)
    filepath = pm.outdir()+'/CaMa_out_'+exp_name+dahour_str+'/'
    base_date = datetime(1871,1,1,0)
    sta_step = int((stime - base_date).days*24 - 1)
    list_time = [sta_step + t for t in range(time_len)]


    # dahour_str = '{:02d}'.format(pm.dahour())
    # filepath = pm.outdir()+'/'+pm.CaMa_out()+pm.dahour()+'/'
    
    def read_bin(ctime,yyyy,enum,chr_type):
        ens = '{:03d}'.format(enum)
        fname = filepath+ctime+chr_type+ens+'/'+var+yyyy+'.bin'
        with open(fname, 'rb') as file:
            data_bin = np.fromfile(file, dtype=np.float32)
            data_re  = np.reshape(data_bin,(dahour,1320,1500))
        return data_re
    
    nhour = 0 
    current_time = stime
    while current_time < etime:
        ctime = current_time.strftime('%Y%m%d%H')
        nyyyy = current_time.strftime('%Y')
        for enum in range(1,1+mean_num):
            data_file[nhour:nhour+dahour,:,:,int(enum-1)] = read_bin(ctime,nyyyy,enum,type_file)
        current_time += time_step
        nhour = nhour+dahour
    
    data_min = np.nanmin(data_file,axis=3)
    data_max = np.nanmax(data_file,axis=3)
    data_mean = np.nanmean(data_file,axis=3)
    
    def make_netcdf_file(name,time_list):
        ncfile = nc.Dataset(savepath+name+".nc", 'w', format='NETCDF4')
        ncfile.createDimension('time', None)
        ncfile.createDimension('lat', 1320) ; ncfile.createDimension('lon', 1500)
        nctime = ncfile.createVariable('time', np.dtype('double').char, ('time'))
        nctime.standard_name = 'time' ; nctime.long_name = 'Time'
        nctime.units = 'hours since 1871-01-01 00:00:00' ; nctime.calendar = 'proleptic_gregorian'
        nctime.axis = 'T'
        nctime[:] = np.array(time_list)
        nclat = ncfile.createVariable('lat', np.dtype('float32').char, ('lat'))
        nclat.standard_name = 'latitude' ; nclat.long_name = 'Latitude'
        nclat.units = 'degrees_north' ; nclat.axis = 'Y'
        nclat[:] = np.linspace(46, 24, 1320)
        nclon = ncfile.createVariable('lon', np.dtype('float32').char, ('lon'))
        nclon.standard_name = 'longitude' ; nclon.long_name = 'Longitude'
        nclon.units = 'degrees_east' ; nclon.axis = 'X'
        nclon[:] = np.linspace(123, 148, 1500)
        return ncfile
    
    def write_netcdf(out_data, ncfile, var, lname, sname, vunit, vfill):
        ncvar = ncfile.createVariable(var, np.dtype('float32').char, ('time', 'lat', 'lon'))
        ncvar.standard_name = sname ; ncvar.long_name = lname ; ncvar.units = vunit
        ncvar._fillvalu = vfill
        ncvar[:,:,:] = out_data
        ncfile.close

    filename = '01_20'

    cdf_file = make_netcdf_file('data'+type_file+'_min'+filename,list_time)
    if var == 'outflw':
        write_netcdf(data_min,cdf_file,'Outflw','River outflow','river outflow','kg m-2 s-1', 1.e+20)
    else:
        write_netcdf(data_min,cdf_file,'Rivdph','River depth','river depth','m', 1.e+20)

    cdf_file = make_netcdf_file('data'+type_file+'_max'+filename,list_time)
    if var == 'outflw':
        write_netcdf(data_max,cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    else:
        write_netcdf(data_max,cdf_file,'Rivdph','River depth','river depth','m', 1.e+20)

    cdf_file = make_netcdf_file('data'+type_file+'_mean'+filename,list_time)
    if var == 'outflw':
        write_netcdf(data_mean,cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    else:
        write_netcdf(data_mean,cdf_file,'Rivdph','River depth','river depth','m', 1.e+20)

generate_nc('rivdph','A')
generate_nc('rivdph','C')
generate_nc('outflw','A')
generate_nc('outflw','C')

















