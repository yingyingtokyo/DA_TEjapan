import numpy as np
import netCDF4 as nc

def read_bin(dahour):
    fname = '/data42/yingying/HydroDA/out_valid/CaMa_out_ils01/2019101210A010/outflw2019.bin'
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(dahour,1320,1500))
    return data_re
data = read_bin(1)
print(data[:,639,590])
print(data[:,623,804])


def read_ens_nc(file):
    nf = nc.Dataset(file,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    print(varname)
    lat = np.array(nf.variables[varname[0]][:])
    lon = np.array(nf.variables[varname[1]][:])
    var = np.array(nf.variables[varname[3]][:])
    var = np.where(var<-100,np.nan,var)
    return lat,lon,var
output_path = '/data42/yingying/HydroDA/out_valid/CaMa_out_ils01/2019101210A010/'
lat,lon,dis = read_ens_nc(output_path+'o_outflw2019.nc')
print(dis[:,639,590])

