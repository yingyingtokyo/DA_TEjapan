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

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,BoundaryNorm
import matplotlib.cm as cm
import imageio

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.io import shapereader

#%%
def generate_file(var,exp_name):
# var  = 'rivdph'
# var  = 'outflw'
mean_num = 10
# exp_name = "tej"
exp_name = "ils"
savepath = '/work/a06/yingying/Code/plot/exp_'+exp_name+'/'+var+'/'
os.makedirs(savepath, exist_ok=True)

# 100100-101823
day1 = 1
day2 = 18
mon  = 1
time_len = (day2-day1+1)*24
ind = 0

da_path = '/work/a06/yingying/da_dph2019/exp_'+exp_name+'/assim_out/'
def read_da(day,hr):
    day_n    = '{:02d}'.format(day)
    hr_n     = '{:02d}'.format(hr)
    filename = day_n+hr_n+'_xa.bin'
    with open(da_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
        day_array = day_array.reshape(1320,1500)
    f.close()
    return day_array

off_data = np.full((time_len,1320,1500),np.nan)
for day in range(day1,day2+1):
    for hr in range(1,25):
        off_data[ind,:,:] = read_da(day,hr)
        off_data[ind,:,:] = np.where(off_data[ind,:,:]==0,np.nan,off_data[ind,:,:])
        ind = ind+1
print(off_data[:,583,1005])

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


cdf_file = make_netcdf_file('dataF_mean01_20')
write_netcdf(off_data[:,:,:],cdf_file,'Rivdph','Rivdph outflow','river outflow','m', 1.e+20)


