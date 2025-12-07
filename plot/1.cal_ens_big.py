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

import sys
params_path = '/data42/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params as pm

#%%
def generate_nc(dahour,var,ens_start,exp_name,type_file):
    dahour_str = '{:02d}'.format(dahour) 
    mean_num = 20
    savepath = '/work/a06/yingying/Code/plot/exp_'+exp_name+dahour_str+'/'+var+'/'
    os.makedirs(savepath, exist_ok=True)
    
    # 100100-101823
    day1 = 1
    day2 = 18
    mon  = 1
    time_len = (day2-day1+1)*24
    data_file = np.full((time_len,1320,1500,mean_num),np.nan)
    filepath = pm.outdir()+'/CaMa_out_'+exp_name+dahour_str+'/'
    # dahour_str = '{:02d}'.format(pm.dahour())
    # filepath = pm.outdir()+'/'+pm.CaMa_out()+pm.dahour()+'/'
    
    def read_bin(mon,date,hour,enum,chr_type):
        month = "{:02d}".format(mon)
        day   = "{:02d}".format(date)
        hr    = "{:02d}".format(hour)
        ens   = "{:03d}".format(enum)
        fname = filepath+'2019'+month+day+hr+chr_type+ens+'/'+var+'2019.bin'
        with open(fname, 'rb') as file:
            data_bin = np.fromfile(file, dtype=np.float32)
            data_re  = np.reshape(data_bin,(dahour,1320,1500))
        return data_re
    
    hour = 0 
    for date in range(0,int(time_len/dahour)):
        start_day = date*dahour
        end_day   = start_day+dahour
        day = int(date//int(24/dahour))+1
        for enum in range(ens_start,ens_start+mean_num):
            ens = enum-ens_start
            data_file[start_day:end_day,:,:,ens] = read_bin(mon,day,hour,enum,type_file)
        hour = hour+dahour
        if hour>23:
            hour = hour-24
    
    data_min = np.nanmin(data_file,axis=3)
    data_max = np.nanmax(data_file,axis=3)
    data_mean = np.nanmean(data_file,axis=3)
    
    def make_netcdf_file(name):
        ncfile = nc.Dataset(savepath+name+".nc", 'w', format='NETCDF4')
        ncfile.createDimension('time', None)
        ncfile.createDimension('lat', 1320) ; ncfile.createDimension('lon', 1500)
        nctime = ncfile.createVariable('time', np.dtype('double').char, ('time'))
        nctime.standard_name = 'time' ; nctime.long_name = 'Time'
        nctime.units = 'hours since 2019-10-01 00:00:00' ; nctime.calendar = 'proleptic_gregorian'
        nctime.axis = 'T'
        nctime[:] = np.arange(0,time_len,1)
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

    if ens_start == 1:
        filename = '01_10'
    else:
        filename = '10_20'

    cdf_file = make_netcdf_file('data'+type_file+'_min'+filename)
    if var == 'outflw':
        write_netcdf(data_min[:,:,:],cdf_file,'Outflw','River outflow','river outflow','kg m-2 s-1', 1.e+20)
    else:
        write_netcdf(data_min[:,:,:],cdf_file,'Rivdph','River depth','river depth','m', 1.e+20)

    cdf_file = make_netcdf_file('data'+type_file+'_max'+filename)
    if var == 'outflw':
        write_netcdf(data_max[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    else:
        write_netcdf(data_max[:,:,:],cdf_file,'Rivdph','River depth','river depth','m', 1.e+20)

    cdf_file = make_netcdf_file('data'+type_file+'_mean'+filename)
    if var == 'outflw':
        write_netcdf(data_mean[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    else:
        write_netcdf(data_mean[:,:,:],cdf_file,'Rivdph','River depth','river depth','m', 1.e+20)

#generate_nc(1,'outflw',1,'ils','A')
#generate_nc(1,'outflw',11,'ils','A')
#generate_nc(1,'outflw',1,'ils','C')
#generate_nc(1,'outflw',11,'ils','C')
#generate_nc(24,'outflw',11,'ils','A')

generate_nc(1,'rivdph',1,'ils','A')
generate_nc(1,'rivdph',1,'ils','C')

















