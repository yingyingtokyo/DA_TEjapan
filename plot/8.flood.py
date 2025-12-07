#%%
import numpy as np
import netCDF4 as nc
import os
import re
import pandas as pd
import glob
import rasterio
from datetime import datetime, timedelta
import pandas as pd
import math

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,BoundaryNorm
import matplotlib.cm as cm
from matplotlib import colors
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
# import imageio
# import pygeodesy
# from pygeodesy.ellipsoidalKarney import LatLon
import folium
# from shapely.geometry import Point,Polygon
# from folium.plugins import FastMarkerCluster

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.io import shapereader
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from scipy.signal import correlate
import matplotlib.colors as mcolors

import sys
params_path = '/data42/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_ils01 as pm

camadir  = pm.CaMa_dir()
inputdir = pm.plot_dir()
expname  = pm.runname(pm.mode())
dahour   = pm.dahour()
dahour_str = '{:02d}'.format(dahour)
exp_plotdir = inputdir + '/exp_' + expname + dahour_str + '/'
obsdir = pm.obs_dir(pm.real_mode())

pm_path = pm.__file__
# simulation start from 20191007 but analysis start from 20191010
start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
syyyy='%04d' % (start_year)
if syyyy == '2019':
    tstart = int(24*3/dahour)
else:
    tstart = 0
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

path_alloc   = pm.alloc_dir()
name_alloc   = '/wlv_2019h.xlsx'
df     = pd.read_excel(path_alloc+name_alloc)
df.rename(columns={df.columns[1]: 'latitude', df.columns[2]: 'longitude'}, inplace=True)
id_info =  np.array(df.iloc[:,0])
lat_obs =  np.array(df.iloc[:,1])
lon_obs =  np.array(df.iloc[:,2])
ix     = np.array(df.iloc[:,8])-1  # lon_ind
iy     = np.array(df.iloc[:,9])-1  # lat_ind

def read_bin2d(fname):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(1320,1500))
    return data_re
def read_bin3d(fname,layer):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(layer,1320,1500))
    return data_re
def read_bin_int(fname,layer):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.int32)
        data_re  = np.reshape(data_bin,(layer,1320,1500))
    return data_re
def read_bin_station(fname):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
    return data_bin
def read_nc(fname,layer):
    nf = nc.Dataset(fname,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    lat = np.array(nf.variables[varname[1]][:])
    lon = np.array(nf.variables[varname[2]][:])
    var = np.array(nf.variables[varname[layer]][:])
    var = np.array(np.where(var>10**7,np.nan,var))
    return var   

# from 2019100100
def read_prep(file):
    nf = nc.Dataset(file,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    print(varname)
    lat = np.array(nf.variables[varname[1]][:])
    lon = np.array(nf.variables[varname[2]][:])
    var = np.array(nf.variables[varname[3]][:])
    var = var*3600 # unit: mm/h 
    var = np.array(np.where(var<-900,np.nan,var))
    return lat,lon,var
lat,lon,rainf = read_prep('/data42/yingying/data/matsiro_frc/tej_rain2019.nc')

lat_min,lat_max,lon_min,lon_max = 34,37,137,141
#%% read flood extend data
# read sim/DA output
def read_flood(dhour):
    da_data = []
    sim_data = []
    SAR_time = datetime(2019,10,12,21)
    time_step = timedelta(hours=dhour)
    ctime = SAR_time - time_step
    cyyyy='%04d' % (ctime.year)
    cmm='%02d' % (ctime.month)
    cdd='%02d' % (ctime.day)
    chh='%02d' % (ctime.hour)
    ctime = cyyyy + cmm + cdd + chh
    chour = math.floor(21/dahour) 
    diff_hour = 21%dahour    # SAR: 2019101221 
    for ens in range(0,pm.ens_mem()):
        ens_str = '{:03d}'.format(ens+1)
        if os.path.exists(pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'A'+ens_str+'/o_fldfrc'+syyyy+'.nc'):
            fname_da  = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'A'+ens_str+'/o_fldfrc'+syyyy+'.nc'
            fname_sim = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'C'+ens_str+'/o_fldfrc'+syyyy+'.nc'
        else:
            chh = '{:02d}'.format(chour)
            ctime = cyyyy + cmm + cdd + chh
            fname_da  = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'A'+ens_str+'/o_fldfrc'+syyyy+'.nc'
            fname_sim = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'C'+ens_str+'/o_fldfrc'+syyyy+'.nc'
        print(fname_da)
        flood_da  = read_nc(fname_da,3)
        flood_sim = read_nc(fname_sim,3)
        lon_region_ind = np.where((lon>=lon_min)&(lon<=lon_max))[0]
        lon_region = lon[lon_region_ind]
        lat_region_ind = np.where((lat>=lat_min)&(lat<=lat_max))[0]
        lat_region = lat[lat_region_ind]
        flood_da_region  = flood_da[dhour,lat_region_ind[:,None],lon_region_ind]
        flood_sim_region = flood_sim[dhour,lat_region_ind[:,None],lon_region_ind]
        da_data.append(flood_da_region)
        sim_data.append(flood_sim_region)
    da_all = np.nanmean(np.stack(da_data,axis=0),axis=0)
    sim_all = np.nanmean(np.stack(sim_data,axis=0),axis=0)
    return da_all,sim_all,lon_region,lat_region,diff_hour
    
# [lat,lon,time_ahead]
flood_da_test,flood_sim_test,lon_region,lat_region,diff_hour = read_flood(21)
num = math.floor(23/dahour)
flood_da = np.full((num,len(lat_region),len(lon_region)),np.nan)
flood_sim = np.full((num,len(lat_region),len(lon_region)),np.nan)
for ind in range(0,num):
    if diff_hour == 0: 
        dhour=dahour+ind*dahour
    else:
        dhour=diff_hour+ind*dahour
    flood_da_each,flood_sim_each,lon_region,lat_region,diff_hour = read_flood(dhour)
    flood_da[ind,:,:] = flood_da_each
    flood_sim[ind,:,:] = flood_sim_each
    
filepath = inputdir+'/exp_'+expname+dahour_str+'/flood/'
os.makedirs(filepath,exist_ok=True)
np.save(filepath+'flood_da.npy',flood_da)
np.save(filepath+'flood_sim.npy',flood_sim)

