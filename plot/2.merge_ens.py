#%%
import numpy as np
import netCDF4 as nc
import os
import re
def generate_nc(dahour,var,type_int,exp_name):
    dahour_str = '{:02d}'.format(dahour)
    filepath = "/work/a06/yingying/Code/plot/"+exp_name+dahour_str+"/"+var+'/'
    mean_list = ['data'+type_int+'_mean01_10','data'+type_int+'_mean10_20']
    min_list = ['data'+type_int+'_min01_10','data'+type_int+'_min10_20']
    max_list = ['data'+type_int+'_max01_10','data'+type_int+'_max10_20']
    day2     = 18
    day1     = 1
    time_len = (day2-day1+1)*24
    min_ens  = np.zeros((time_len,1320,1500,len(min_list)))
    max_ens  = np.zeros((time_len,1320,1500,len(max_list)))
    mean_ens  = np.zeros((time_len,1320,1500,len(mean_list)))
    
    def read_ens_nc(file):
        nf = nc.Dataset(file,'r')
        varname = nf.variables.keys()
        varname = list(varname)
        print(varname)
        lat = np.array(nf.variables[varname[0]][:])
        lon = np.array(nf.variables[varname[1]][:])
        var = np.array(nf.variables[varname[3]][:])
        var = np.where(var<-100,np.nan,var)
        var = np.where(var>10**8,np.nan,var)
        return lat,lon,var
    
    ind=0
    for min_n in min_list:
        lat,lon,min_ens[:,:,:,ind] = read_ens_nc(filepath+min_n+'.nc')
        ind = ind+1
    min_data = np.nanmin(min_ens,axis=3)
    
    ind=0
    for max_n in max_list:
        lat,lon,max_ens[:,:,:,ind] = read_ens_nc(filepath+max_n+'.nc')
        ind = ind+1
    max_data = np.nanmax(max_ens,axis=3)
    
    ind=0
    for mean_n in mean_list:
        lat,lon,mean_ens[:,:,:,ind] = read_ens_nc(filepath+mean_n+'.nc')
        ind = ind+1
    mean_data = np.nanmean(mean_ens,axis=3)
    
    def make_netcdf_file(name):
        ncfile = nc.Dataset(filepath+name+".nc", 'w', format='NETCDF4')
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
    
    cdf_file = make_netcdf_file('data'+type_int+'_min01_20')
    write_netcdf(min_data[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    
    cdf_file = make_netcdf_file('data'+type_int+'_max01_20')
    write_netcdf(max_data[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)
    
    cdf_file = make_netcdf_file('data'+type_int+'_mean01_20')
    write_netcdf(mean_data[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)

# dahour = 3
# var = 'rivdph'
# var = 'outflw'
# type_int = 'C'
#type_int = 'A'
#exp_name = 'exp_wid'
#exp_name = 'exp_ils'
#exp_name = 'exp_hgt'
# exp_name = 'exp_man'
# generate_nc(3,'rivdph','C','exp_man')
# generate_nc(3,'rivdph','A','exp_man')
# generate_nc(3,'outflw','C','exp_man')
# generate_nc(3,'outflw','A','exp_man')
generate_nc(3,'rivdph','C','exp_man')
generate_nc(3,'rivdph','A','exp_man')
generate_nc(3,'outflw','C','exp_man')
generate_nc(3,'outflw','A','exp_man')
#generate_nc(24,'outflw','C','exp_man')
#generate_nc(24,'outflw','A','exp_man')
# generate_nc(24,'rivdph','C','exp_ils')
# generate_nc(24,'rivdph','A','exp_ils')
# generate_nc(24,'outflw','C','exp_ils')
# generate_nc(24,'outflw','A','exp_ils')
#generate_nc(12,'rivdph','C','exp_ils')
#generate_nc(12,'rivdph','A','exp_ils')
#generate_nc(12,'outflw','C','exp_ils')
#generate_nc(12,'outflw','A','exp_ils')







