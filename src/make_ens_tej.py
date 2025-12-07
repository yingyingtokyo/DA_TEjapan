# -*- coding: utf-8 -*-
"""
Make TE-japan as the forcing dataset of CAMA-Flood
@author: LIU YINGYING
"""

#var_name = "elevtn"
#var_name = "rivwth_gwdlr"
#var_name = "rivhgt"
#var_name = "rivman"
var_name = "runoff"
days=30

#%%
import os
import numpy as np
import netCDF4 as nc
import numpy.random as rd

#%% read nc
ens_num  = 10
savepath = '/work/a06/yingying/CaMa_v411/cmf_v411_pkg/inp/tej'
dis_data = np.full((720,1320,1500,ens_num),0,dtype=np.float32)

#%% read nc
def read_nc():
    file = '/work/a06/yingying/Runoff/201910/TE_tej_MSM.nc'
    # check if the existence of data
    if os.path.exists(file):
        nf = nc.Dataset(file,'r')
        varname = list(nf.variables.keys())
        var = np.array(nf.variables[varname[3]][:,:,:])
        # var = var*1000/3600  #ERA5 # kg/m^2/s -> mm/day
        var = var*86400  #tej
        var = np.where(var<-100,0,var)
        var = np.where(var>10**8,0,var)
    return var
data_re = read_nc()

def make_rand(): #used
    # prepare random numbers for ensemable simulation
    for lat_num in range(0,1320):
        if np.all(data_re[:,lat_num,:]<10**(-12))==True:
            continue
        for lon_num in range(0,1500):
            if np.all(data_re[:,lat_num,lon_num]<10**(-12))==True:
                continue
            else:
                if var_name == "runoff":
                    rand = np.random.lognormal(mean=1,sigma=1.2,size=ens_num)  # runoff
                    for num in range(0,ens_num):
                        dis_data[:,lat_num,lon_num,num] = rand[num]*data_re[:,lat_num,lon_num]
    return dis_data 

runoff = make_rand()

for ens in range(0,ens_num):
    str_ens = '{:03d}'.format(ens+1)
    var_data = runoff[:,:,:,ens]
    #var_data = var[:,:,:]
    for mon in range(0,12):
        str_mon = '{:02d}'.format(mon+1)
        for ind in range(0,days):
            var_each = np.array(var_data[ind*24:(ind+1)*24,:,:], dtype=np.float32)
            var_bin = var_each.tobytes()
            # write .bin
            fileday ='{:02d}'.format(ind+1)
            filebin  = savepath+'/Runoff_2019'+str_mon+fileday+str_ens+'.bin'
            with open(filebin,'wb') as file:
                file.write(var_bin)
    os.rename(savepath+'/Runoff_20190229'+str_ens+'.bin',savepath+'/Runoff_20190131'+str_ens+'.bin')
    os.rename(savepath+'/Runoff_20190230'+str_ens+'.bin',savepath+'/Runoff_20190331'+str_ens+'.bin')
    shutil.copy(savepath+'/Runoff_20190601'+str_ens+'.bin',savepath+'/Runoff_20190531'+str_ens+'.bin')
    shutil.copy(savepath+'/Runoff_20190801'+str_ens+'.bin',savepath+'/Runoff_20190731'+str_ens+'.bin')
    shutil.copy(savepath+'/Runoff_20190901'+str_ens+'.bin',savepath+'/Runoff_20190831'+str_ens+'.bin')
    shutil.copy(savepath+'/Runoff_20191101'+str_ens+'.bin',savepath+'/Runoff_20191031'+str_ens+'.bin')
    shutil.copy(savepath+'/Runoff_20190101'+str_ens+'.bin',savepath+'/Runoff_20191231'+str_ens+'.bin')




