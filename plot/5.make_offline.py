import numpy as np
import netCDF4 as nc
import os

#exp  = 'exp_tej' # change runoff
#exp  = 'exp_ils' # change rainfall
#exp  = 'exp_wid' # change wid
#exp  = 'exp_hgt' # change hgt
exp  = 'exp_man' # change manning

dahour = 3
dahour_str = '{:02d}'.format(dahour)
day1 = 1
day2 = 18
time_len = int((day2-day1+1)*24/dahour)
ens_num  = 20

var_name = 'rivdph' 
#var_name = 'outflw' 
filepath = "/work/a06/yingying/Code/plot/"+exp+dahour_str+"/"+var_name+'/'

def read_bin(map_path,filename):
    with open(map_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
        day_array = day_array.reshape(1320,1500)
        day_array = np.where(day_array<-999,np.nan,day_array)
    return day_array

# offline DA
offpath = '/work/a06/yingying/camada/HydroDA/src/assim_'+exp[-3:]+dahour_str+'/ens_xa/meanA/201901'
dataF = np.full((int((day2-day1+1)*24/dahour),1320,1500),np.nan)
shour = 0
for day in range(0,time_len):
    date  = int(day//int(24/dahour))+1
    fday = '{:02d}'.format(date)
    fhour = '{:02d}'.format(shour)
    dataF[day,:,:]=read_bin(offpath,fday+fhour+'_xa.bin')
    shour = shour + dahour
    if shour>23:
        shour = shour-24
print(dataF[:,612,974])

def make_netcdf_file(name):
    ncfile = nc.Dataset(filepath+name+".nc", 'w', format='NETCDF4')
    ncfile.createDimension('lat', 1320) ; ncfile.createDimension('lon', 1500)
    ncfile.createDimension('time', None)
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


cdf_file = make_netcdf_file('dataF')
write_netcdf(dataF[:,:,:],cdf_file,'Outflw','Rivdph outflow','river outflow','kg m-2 s-1', 1.e+20)

