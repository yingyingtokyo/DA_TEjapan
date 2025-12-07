mport numpy as np
import netCDF4 as nc
import os
import re
import pandas as pd
import glob
import rasterio
from datetime import datetime, timedelta
import pandas as pd

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
from matplotlib.patches import Rectangle

camadir  = pm.CaMa_dir()
inputdir = pm.plot_dir()
expname  = pm.runname(pm.mode())
dahour   = pm.dahour()
dahour_str = '{:02d}'.format(dahour)
exp_plotdir = inputdir + '/exp_' + expname + dahour_str + '/'
obsdir = pm.obs_dir(pm.real_mode())

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
rainf = rainf[6*24:16*24,:,:]

lat_min,lat_max,lon_min,lon_max = 30,46,128,149
def read_plot_tif(filename,lat_min,lat_max,lon_min,lon_max):
    up_area  = np.zeros(((lat_max-lat_min)*3600,(lon_max-lon_min)*3600)) 
    lat_all  = np.zeros((lat_max-lat_min)*3600) 
    lon_all  = np.zeros((lon_max-lon_min)*3600)
    filepath = pm.DA_dir()+'/src/plot/file/map/'+filename
    tiff_files = glob.glob(os.path.join(filepath, '*.tif'))
    for filenum in range(0,len(tiff_files)):
        filen = tiff_files[filenum]
        file_lon = int(filen[51:54])
        file_lat = int(filen[48:50])
        if (file_lat-lat_min>=0) &(file_lon-lon_min>=0) &(file_lat-lat_max<0) &(file_lon-lon_max<0):           
            with rasterio.open(filen) as dataset:
                tiff_data = dataset.read(1)
                # be care !! when drawing the figure should convert latitude tiff[::-1,:]
                up_area[(file_lat-lat_min)*3600:(file_lat-lat_min+1)*3600,(file_lon-lon_min)*3600:(file_lon-lon_min+1)*3600] = tiff_data[::-1,:]
    return up_area
up_plot = read_plot_tif('upa',30,46,128,149)
#%% lower resolution
size=120
uparea_filter = np.where(up_plot<10**2,np.nan,up_plot)
def low_resolution(var,size):
    var_low = np.full((int((lat_max-lat_min)*3600/size),int((lon_max-lon_min)*3600/size)),np.nan)
    for row in range(0,np.shape(var_low)[0]):
        if (np.all(np.isnan(var[row*size:(row+1)*size,:]))==1):
            continue
        for col in range(0,np.shape(var_low)[1]):
            var_low[row,col] = np.nanmean(var[row*size:(row+1)*size,col*size:(col+1)*size])
    return var_low    
up_low = low_resolution(uparea_filter,size) # for figure ploting
lon_map_low = np.linspace(lon_min,lon_max,int((lon_max-lon_min)*3600/size))
lat_map_low = np.linspace(lat_min,lat_max,int((lat_max-lat_min)*3600/size))
lon2d_low_map, lat2d_low_map = np.meshgrid(lon_map_low, lat_map_low)

# #%% read small tiff data
lat_min,lat_max,lon_min,lon_max = 34,37,137,141
up_area = read_plot_tif('upa',34,37,137,141)
#%% lower resolution
size=10
def low_resolution(var,size):
    var_low = np.full((int((lat_max-lat_min)*3600/size),int((lon_max-lon_min)*3600/size)),np.nan)
    for row in range(0,np.shape(var_low)[0]):
        if (np.all(np.isnan(var[row*size:(row+1)*size,:]))==1):
            continue
        for col in range(0,np.shape(var_low)[1]):
            var_low[row,col] = np.nanmean(var[row*size:(row+1)*size,col*size:(col+1)*size])
    return var_low    
uparea_low = low_resolution(up_area,size)
uparea_low = np.where(uparea_low<10,np.nan,uparea_low)
lon_tiff = np.linspace(lon_min,lon_max,int((lon_max-lon_min)*3600/size))
lat_tiff = np.linspace(lat_min,lat_max,int((lat_max-lat_min)*3600/size))


#%% plot flood extend
filepath = '/data50/yingying/flood/'
def read_cov_tif(filename):
    filen = filepath + filename
    with rasterio.open(filen) as dataset:
        tiff_data = dataset.read(1)
        # tiff_data = tiff_data[::-1,:]
        tiff_data = np.where(tiff_data==0,np.nan,tiff_data)
        col = dataset.width      
        row = dataset.height     
        bounds = dataset.bounds
        min_lon, min_lat, max_lon, max_lat = bounds
        print(f"Longitude Range: {min_lon} to {max_lon}")
        print(f"Latitude Range: {min_lat} to {max_lat}")
        print('--------')
    return col,row,min_lon,max_lon,min_lat,max_lat,tiff_data
# 20191012T2042 UST
col1,row1,min_lon1,max_lon1,min_lat1,max_lat1,fld_reg1 = read_cov_tif('area_1_result.tif')  # 138-138.4, 36.6-36.8
col2,row2,min_lon2,max_lon2,min_lat2,max_lat2,fld_reg2 = read_cov_tif('area_2_result.tif')  # 139.3-139.8, 34.45-35.7

# col1,row1,min_lon1,max_lon1,min_lat1,max_lat1,fld_reg1 = read_cov_tif('binary_output_without_water_area2.tif')  # 137.4-140.6, 35.9-37.8
# col2,row2,min_lon2,max_lon2,min_lat2,max_lat2,fld_reg2 = read_cov_tif('binary_output_without_water_coast.tif')  # 138.1-141.2, 34.5-36.4

lon_map1 = np.linspace(min_lon1,max_lon1,col1)
lat_map1 = np.linspace(min_lat1,max_lat1,row1)
lon2d_map1, lat2d_map1 = np.meshgrid(lon_map1, lat_map1)
lon_map2 = np.linspace(min_lon2,max_lon2,col2)
lat_map2 = np.linspace(min_lat2,max_lat2,row2)
lon2d_map2, lat2d_map2 = np.meshgrid(lon_map2, lat_map2)

d flood extend data
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
    for ens in range(0,pm.ens_mem()):
        ens_str = '{:03d}'.format(ens+1)
        fname_da  = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'A'+ens_str+'/o_fldfrc'+syyyy+'.nc'
        fname_sim = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'+ctime+'C'+ens_str+'/o_fldfrc'+syyyy+'.nc'
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
        break
    # da_all = np.nanmean(np.stack(da_data,axis=0),axis=0)
    # sim_all = np.nanmean(np.stack(sim_data,axis=0),axis=0)
    # print(np.shape(da_all)) # [180,240]
    return lon_region,lat_region
    
# [lat,lon,time_ahead]
lon_region,lat_region = read_flood(21)
# flood_da = np.full((len(lat_region),len(lon_region),23),np.nan)
# flood_sim = np.full((len(lat_region),len(lon_region),23),np.nan)
# for dhour in range(1,24):
#     flood_da_each,flood_sim_each,lon_region,lat_region = read_flood(dhour)
#     flood_da[:,:,dhour-1] = flood_da_each
#     flood_sim[:,:,dhour-1] = flood_sim_each
    

filepath = inputdir+'/exp_'+expname+dahour_str+'/'
# np.save(filepath+'flood'+'/flood_da.npy',flood_da)
# np.save(filepath+'flood'+'/flood_sim.npy',flood_sim)
flood_da = np.load(filepath+'flood'+'/flood_da.npy')
flood_sim = np.load(filepath+'flood'+'/flood_sim.npy')

def find_sim_region(min_lon,max_lon,min_lat,max_lat):
    lat_ind = np.where((lat_region>=min_lat)&(lat_region<=max_lat))[0]
    lon_ind = np.where((lon_region>=min_lon)&(lon_region<=max_lon))[0]
    lat_reg = lat_region[np.nanmin(lat_ind):np.nanmax(lat_ind)+1]
    lon_reg = lon_region[np.nanmin(lon_ind):np.nanmax(lon_ind)+1]
    lon2d_map, lat2d_map = np.meshgrid(lon_reg, lat_reg)
    sim_fld_reg = flood_sim[np.nanmin(lat_ind):np.nanmax(lat_ind)+1,np.nanmin(lon_ind):np.nanmax(lon_ind)+1,:]
    da_fld_reg  = flood_da[np.nanmin(lat_ind):np.nanmax(lat_ind)+1,np.nanmin(lon_ind):np.nanmax(lon_ind)+1,:]
    return lon2d_map, lat2d_map,sim_fld_reg,da_fld_reg
lon2d_sim1,lat2d_sim1,sim_fld_reg1,da_fld_reg1 = find_sim_region(min_lon1,max_lon1,min_lat1,max_lat1)
lon2d_sim2,lat2d_sim2,sim_fld_reg2,da_fld_reg2 = find_sim_region(min_lon2,max_lon2,min_lat2,max_lat2)

# upscaling for SAR data
window_size = 46   #93 grids, depends on resolution
def upscale_sat(fld_reg,sim_fld_reg,da_fld_reg):
    fld_upscale = np.full((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1]),np.nan)
    for i in range(0,np.shape(sim_fld_reg)[0]): # lat
        for j in range(0,np.shape(sim_fld_reg)[1]): # lon
            sim = sim_fld_reg[i,j]
            da  = da_fld_reg[i,j]
            lat_i1 = i*window_size*2-window_size
            lat_i2 = i*window_size*2+window_size
            lon_i1 = j*window_size*2-window_size
            lon_i2 = j*window_size*2+window_size
            if lat_i1<0:
                lat_i1 = 0
            if lon_i1<0:
                lon_i1 = 0            
            if lat_i2>np.shape(fld_reg)[0]:
                lat_i2 = np.shape(fld_reg)[0]    
            if lon_i2>np.shape(fld_reg)[1]:
                lon_i2 = np.shape(fld_reg)[1]
            # print(lat_i1,lat_i2,lon_i1,lon_i2]))
            fld_data = fld_reg[lat_i1:lat_i2,lon_i1:lon_i2]
            fld_upscale[i,j] = np.nansum(fld_data)/((lat_i2-lat_i1)*(lon_i2-lon_i1))
    return fld_upscale
fld_up1 = upscale_sat(fld_reg1,sim_fld_reg1,da_fld_reg1)
fld_up2 = upscale_sat(fld_reg2,sim_fld_reg2,da_fld_reg2)

# region compare
def cal_comp_fld(fld_up,sim_fld_reg,da_fld_reg):
    cal_sim_fld = np.zeros((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1]))
    cal_da_fld  = np.zeros((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1]))
    num_fld     = np.full((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1]),np.nan)
    cal_fld     = np.full((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1]),np.nan)
    size_fld    = np.full((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1]),np.nan)
    for i in range(0,np.shape(sim_fld_reg)[0]): # lat
        for j in range(0,np.shape(sim_fld_reg)[1]): # lon
            # simulation condition
            if (fld_up[i,j]<10**(-8)) & (sim_fld_reg[i,j]>10**(-8)):  # only sim
                cal_sim_fld[i,j] = 1
            if (fld_up[i,j]>10**(-8)) & (sim_fld_reg[i,j]>10**(-8)):  # both sim and sat
                cal_sim_fld[i,j] = 2
            if (fld_up[i,j]>10**(-8)) & (sim_fld_reg[i,j]<10**(-8)):  # only sat
                cal_sim_fld[i,j] = 3
            # da condition
            if (fld_up[i,j]<10**(-8)) & (da_fld_reg[i,j]>10**(-8)):  # only da
                cal_da_fld[i,j] = 1
            if (fld_up[i,j]>10**(-8)) & (da_fld_reg[i,j]>10**(-8)):  # both da and sat
               cal_da_fld[i,j] = 2
            if (fld_up[i,j]>10**(-8)) & (da_fld_reg[i,j]<10**(-8)):  # only sat
                cal_da_fld[i,j] = 3
            # 0.5: only sim, 1.5: only da, 2.5: only sat, 3.5: sim+da, 4.5: sim+sat, 5.5: da+sat, 6.5: sim+da+sat, nan: no flood
            if (cal_sim_fld[i,j]==1) & (cal_da_fld[i,j]==0):
                cal_fld[i,j] = 0.5
                size_fld[i,j]= sim_fld_reg[i,j]*100
            if (cal_sim_fld[i,j]==0) & (cal_da_fld[i,j]==1):
                cal_fld[i,j] = 1.5
                size_fld[i,j]= da_fld_reg[i,j]*100
            if (cal_sim_fld[i,j]==3) & (cal_da_fld[i,j]==3):
                cal_fld[i,j] = 2.5                
                size_fld[i,j]= fld_up[i,j]*100
            if (cal_sim_fld[i,j]==1) & (cal_da_fld[i,j]==1):
                cal_fld[i,j] = 3.5            
                size_fld[i,j]= sim_fld_reg[i,j]*100
            if (cal_sim_fld[i,j]==2) & (cal_da_fld[i,j]==3):
                cal_fld[i,j] = 4.5                
                size_fld[i,j]= fld_up[i,j]*100
            if (cal_sim_fld[i,j]==3) & (cal_da_fld[i,j]==2):
                cal_fld[i,j] = 5.5
                size_fld[i,j]= fld_up[i,j]*100
            if (cal_sim_fld[i,j]==2) & (cal_da_fld[i,j]==2):
                cal_fld[i,j] = 6.5      
                size_fld[i,j]= fld_up[i,j]*100
    sim1 = np.shape(np.where(cal_sim_fld==1))[1]
    sim2 = np.shape(np.where(cal_sim_fld==2))[1]
    sim3 = np.shape(np.where(cal_sim_fld==3))[1]
    csi_sim  = sim2/(sim1+sim2+sim3)
    da1 = np.shape(np.where(cal_da_fld==1))[1]
    da2 = np.shape(np.where(cal_da_fld==2))[1]
    da3 = np.shape(np.where(cal_da_fld==3))[1]
    csi_da   = da2/(da1+da2+da3)
    # def cal_csi(cal_sim_fld):
    #     sim1 = np.where(cal_sim_fld==1)
    print('score sim:',sim1,sim2,sim3,csi_sim,'score_da:',da1,da2,da3,csi_da)
    return cal_fld,size_fld,cal_sim_fld,cal_da_fld,csi_sim,csi_da

# from forecast time ahead 1hour to 24hours
def evaluate_flood(fld_up,sim_fld_reg,da_fld_reg):
    cal_sim_fld = np.zeros((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1],np.shape(sim_fld_reg)[2]))
    cal_da_fld  = np.zeros((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1],np.shape(sim_fld_reg)[2]))
    cal_fld     = np.full((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1],np.shape(sim_fld_reg)[2]),np.nan)
    size_fld    = np.full((np.shape(sim_fld_reg)[0],np.shape(sim_fld_reg)[1],np.shape(sim_fld_reg)[2]),np.nan)
    csi_sim = np.full(np.shape(sim_fld_reg)[2],np.nan)
    csi_da  = np.full(np.shape(sim_fld_reg)[2],np.nan)
    for i in range(0,np.shape(sim_fld_reg)[2]):
        cal_fld[:,:,i],size_fld[:,:,i],cal_sim_fld[:,:,i],cal_da_fld[:,:,i],csi_sim[i],csi_da[i] = cal_comp_fld(fld_up,sim_fld_reg[:,:,i],da_fld_reg[:,:,i])
    return cal_fld,size_fld,cal_sim_fld,cal_da_fld,csi_sim,csi_da
cal_fld1,size_fld1,cal_sim_fld1,cal_da_fld1,csi_sim1,csi_da1 = evaluate_flood(fld_up1,sim_fld_reg1,da_fld_reg1)
cal_fld2,size_fld2,cal_sim_fld2,cal_da_fld2,csi_sim2,csi_da2 = evaluate_flood(fld_up2,sim_fld_reg2,da_fld_reg2)

#%% plot flood condition for 1-24 hours ahead
def plot_flood_condition_all():
    fig, axes = plt.subplots(8, 6, figsize=(36, 24), dpi=600,subplot_kw={'projection': ccrs.PlateCarree()})
    axes = axes.flatten()
    def comp_fld(ax1, lon_ret1, lon_ret2,lat_ret1, lat_ret2, lon2d_sim, lat2d_sim, cal_fld, size_fld, plot_idx):
        ax1.set_extent([lon_ret1, lon_ret2,lat_ret1, lat_ret2], ccrs.PlateCarree())
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=1, color='gray', alpha=0, linestyle='-.')
        ax1.coastlines(alpha=1.,linestyle='-',color = "#3F3A3A",lw=0.8,resolution='10m')
        cmap_blue = ListedColormap(['blue']) 
        colors = plt.cm.rainbow(np.linspace(0, 1, 7)) 
        cmap = ListedColormap(colors) 
        bounds = np.arange(0, 8)
        norm = BoundaryNorm(bounds, cmap.N) 
        gl.top_labels   = False
        gl.right_labels = False
        gl.bottom_labels= False
        gl.left_labels = False
        lon_map = np.linspace(lon_min,lon_max,int((lon_max-lon_min)*3600))
        lat_map = np.linspace(lat_min,lat_max,int((lat_max-lat_min)*3600))
        lon_tiff_ind = np.where((lon_map>=lon_ret1)&(lon_map<=lon_ret2))[0]
        lat_tiff_ind = np.where((lat_map>=lat_ret1)&(lat_map<=lat_ret2))[0]
        lat_tiff_reg = lat_map[np.nanmin(lat_tiff_ind):np.nanmax(lat_tiff_ind)+1]
        lon_tiff_reg = lon_map[np.nanmin(lon_tiff_ind):np.nanmax(lon_tiff_ind)+1]
        lon2d_tiff, lat2d_tiff = np.meshgrid(lon_tiff_reg, lat_tiff_reg)
        uparea_reg = up_area[np.nanmin(lat_tiff_ind):np.nanmax(lat_tiff_ind)+1,np.nanmin(lon_tiff_ind):np.nanmax(lon_tiff_ind)+1]
        uparea_reg = np.where(uparea_reg<10,np.nan,0.01)
        ax1.scatter(lon2d_tiff, lat2d_tiff,uparea_reg, c='blue',zorder=3, transform=ccrs.PlateCarree())
        ## ax1.pcolormesh(lon_tiff_reg,lat_tiff_reg,uparea_reg,zorder=2, cmap=cmap_blue, transform=ccrs.PlateCarree())  
        con_plot2 = ax1.scatter(lon2d_sim, lat2d_sim, c=cal_fld, s=size_fld,cmap=cmap, norm=norm, zorder=3, transform=ccrs.PlateCarree())
        if plot_idx == 0:
            def color_bar(l,b,w,h):
              rect = [l,b,w,h]
              cbar_ax = fig.add_axes(rect)
              return cbar_ax
            [l2,b2,w2,h2] = [0.85,0.28,0.03,0.75]
            cbar_ax2 = color_bar(l2,b2,w2,h2)
            cb2 = plt.colorbar(con_plot2, cax=cbar_ax2,orientation="vertical",shrink=0.5)
            cb2.set_ticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
            cb2.set_ticklabels(['only Sim','only DA','only SAR','Sim+DA','Sim+SAR','DA+SAR','Sim+DA+SAR'],fontdict=font_label)
    # 0.5: only sim, 1.5: only da, 2.5: only sat, 3.5: sim+da, 4.5: sim+sat, 5.5: da+sat, 6.5: sim+da+sat, nan: no flood
    for plot_idx in range(0,23):
        comp_fld(axes[plot_idx],lon_ret1, lon_ret2,lat_ret1, lat_ret2, lon2d_sim1, lat2d_sim1, cal_fld1[:,:,plot_idx], size_fld1[:,:,plot_idx], plot_idx)
    axes[24-1].axis('off')
    for plot_idx in range(0,23):
        comp_fld(axes[plot_idx+24],lon_ret3, lon_ret4,lat_ret3, lat_ret4, lon2d_sim2, lat2d_sim2, cal_fld2[:,:,plot_idx], size_fld2[:,:,plot_idx], plot_idx)
    axes[48-1].axis('off')
    plt.tight_layout()
    plt.savefig(exp_plotdir+'6.24hour_flood_condition.jpg', format='jpg',dpi=600)
    plt.close()
plot_flood_condition_all()

