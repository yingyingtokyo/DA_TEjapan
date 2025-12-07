import numpy as np
import netCDF4 as nc

dahour = 24
dahour_str = '{:02d}'.format(dahour)
da_path = '/work/a06/yingying/camada/HydroDA/src/assim_man'+dahour_str+'/ens_xa/assim/'
daname = '2019011300019A_xa.bin'
def read_obs(obs_path,filename):
    with open(obs_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
        day_array = day_array.reshape(1320,1500)
        day_array = np.where(day_array<-900,np.nan,day_array)
    f.close()
    return day_array
data = read_obs(da_path,daname)
print(data[605,993])

def read_bin(map_path,filename):
    with open(map_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
        day_array = day_array.reshape(dahour,1320,1500)
        day_array = np.where(day_array<-999,np.nan,day_array)
    return day_array
sim_path = '/work/a06/yingying/camada/HydroDA/src/CaMa_out_man'+dahour_str+'/2019011300A002/'
#sim = read_bin(sim_path,'outflw2019.bin')
#print(sim[:,605,993])
simmean = np.full((1320,1500,20),np.nan)
for num in range(0,20):
    enum = '{:03d}'.format(num+1)
    simmean[:,:,num] = read_obs('/work/a06/yingying/camada/HydroDA/src/simmean/ils/','C'+enum+'.bin')
simmean_mean = np.nanmean(simmean,axis=2)
print(simmean_mean[603,1034])

obs_path = '/work/a06/yingying/obs/2019UST/wlv/'
#obs = read_obs(obs_path,'2019101212.bin')
#print(obs[605,993])

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
output_path = '/work/a06/yingying/Code/plot/exp_man12/'
lat,lon,disC_mean = read_ens_nc(output_path+'rivdph/dataC_mean01_20.nc')
simmean = np.nanmean(disC_mean,axis=0)
print(simmean[603,1034])
print(disC_mean[:,603,1034])
print(simmean[173,1289])
print(disC_mean[:,173,1289])
print('-----------')
print(simmean[605,993])
print(disC_mean[:,605,993])

exit()
day1 = 1
day2 = 14
time_len  = int((day2-day1+1)*24/dahour)
dataF = np.full((time_len,1320,1500),np.nan)
# offline DA
offpath = '/work/a06/yingying/camada/HydroDA/src/assim_man12/ens_xa/meanA/201901'
shour = 0
for day in range(0,int((day2-day1+1)*24/dahour)):
    date  = int(day//int(24/dahour))+1
    fday = '{:02d}'.format(date)
    fhour = '{:02d}'.format(shour)
    time_each = np.arange(day*dahour,(day+1)*dahour,1)
    dataF[day,:,:]=read_bin(offpath,fday+fhour+'_xa.bin')
    print(time_each[-1],day,fday,fhour)
    print(dataF[day,605,993])
    #print(dataF[day,318,1047])
    print('==============')
    shour = shour + dahour
    if shour>23:
        shour = shour-24

