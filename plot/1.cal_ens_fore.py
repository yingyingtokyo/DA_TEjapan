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


#%%
def generate_nc(dahour,var,exp_name,type_file,leadtime):
    dahour_str = '{:02d}'.format(dahour) 
    leadtime_str = '{:02d}'.format(leadtime)
    mean_num = 20
    savepath = '/work/a06/yingying/Code/fore/exp_'+exp_name+dahour_str+'/'+var+'/'
    os.makedirs(savepath, exist_ok=True)
    
    # 100100-101823
    day1 = 11
    day2 = 14
    mon  = 1
    time_len  = (day2-day1+1)*24
    data_file = np.full((leadtime,1320,1500,mean_num),np.nan)
    filepath  = '/work/a06/yingying/camada/HydroDA/src/CaMa_out_'+exp_name+dahour_str+'/'

    def read_nc(mon,date,hour,enum,chr_type,var):
        month = "{:02d}".format(mon)
        day   = "{:02d}".format(date)
        hr    = "{:02d}".format(hour)
        ens   = "{:03d}".format(enum)
        filenc = filepath+'2019'+month+day+hr+chr_type+ens+'/o_'+var+'2019.nc'
        # read nc file
        nf = nc.Dataset(filenc,'r')
        varname = nf.variables.keys()
        varname = list(varname)
        var = np.array(nf.variables[varname[3]][:])
        var = np.where(var<-1000,np.nan,var)
        var = np.where(var>10**8,np.nan,var)
        return var      
    
    def make_netcdf_file(name):
        ncfile = nc.Dataset(savepath+name+".nc", 'w', format='NETCDF4')
        ncfile.createDimension('time', None)
        ncfile.createDimension('lat', 1320) ; ncfile.createDimension('lon', 1500)
        nctime = ncfile.createVariable('time', np.dtype('double').char, ('time'))
        nctime.standard_name = 'time' ; nctime.long_name = 'Time'
        nctime.units = 'hours since 2019-10-01 00:00:00' ; nctime.calendar = 'proleptic_gregorian'
        nctime.axis = 'T'
        nctime[:] = np.arange(0,24,1)
        nclat = ncfile.createVariable('lat', np.dtype('float32').char, ('lat'))
        nclat.standard_name = 'latitude' ; nclat.long_name = 'Latitude'
        nclat.units = 'degrees_north' ; nclat.axis = 'Y'
        nclat[:] = np.linspace(46, leadtime, 1320)
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

    hour = 0
    for date in range(0,int(time_len/dahour)):
        start_day = date*dahour
        end_day   = start_day+dahour
        day = int(date//int(24/dahour))+day1
        mon_n  = "{:02d}".format(mon)
        date_n = "{:02d}".format(day)
        hour_n = "{:02d}".format(hour)
        
        for enum in range(0,mean_num):
            data_file[:,:,:,enum] = read_nc(mon,day,hour,enum+1,type_file,var)
        data_min = np.nanmin(data_file,axis=3)
        data_max = np.nanmax(data_file,axis=3)
        data_mean = np.nanmean(data_file,axis=3)

        cdf_file = make_netcdf_file('data'+type_file+'_min_'+mon_n+date_n+hour_n)
        write_netcdf(data_min[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    
        cdf_file = make_netcdf_file('data'+type_file+'_max_'+mon_n+date_n+hour_n)
        write_netcdf(data_max[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    
        cdf_file = make_netcdf_file('data'+type_file+'_mean_'+mon_n+date_n+hour_n)
        write_netcdf(data_mean[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
            
        hour = hour+dahour
        if hour>23:
            day  = day+1
            hour = hour-24

generate_nc(1,'rivdph','ils','A',24)
generate_nc(1,'rivdph','ils','A',24)
generate_nc(1,'rivdph','ils','C',24)
generate_nc(1,'rivdph','ils','C',24)
#generate_nc(24,'outflw',11,'ils','A')
# generate_nc(24,'rivdph',1,'man','C')
# generate_nc(24,'rivdph',11,'man','C')

















