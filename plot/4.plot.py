#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:06:49 2024
# draw the river map in each basin

@author: yingyingliu
"""
#%%
import numpy as np
import netCDF4 as nc
import os
import re
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,BoundaryNorm
import matplotlib.cm as cm
import imageio

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.io import shapereader
from sklearn.metrics import mean_squared_error
import matplotlib.colors as mcolors

#%% change the settings
dahour = 6 
dahour_str = '{:02d}'.format(dahour) 
day1 = 1
day2 = 18
time_len = (day2-day1+1)*24
var_dis  = 'outflw' # check discharge
var_dep  = 'rivdph' # check water depth
ens_num  = 20

var_name = var_dis
#var_name = var_dep

if var_name == 'outflw':
    var_ylabel = 'Outflow (kg m-2 s-1)'
    var_dir = 'dis'
else:
    var_ylabel  = 'River Depth (m)'
    var_dir = 'wlv'

#exp  = 'exp_tej' # change runoff
#exp  = 'exp_ils' # change rainfall
#exp  = 'exp_wid' # change wid 
#exp  = 'exp_hgt' # change hgt 
exp  = 'exp_man' # change manning

df  = pd.read_excel("/work/a06/yingying/MLIT/MLIT_FlowGauge_all_info.xlsx",sheet_name='All_locationsV2')
site_name = np.array(df.iloc[:,0])
ref = np.array(df.iloc[:,1])

#%% read the location of observations
file_alloc = 'gauge_good_wlv.txt'
# read allocate information
id_info = list()
ix      = list()
iy      = list()
with open('/work/a06/yingying/obs/alloc/'+file_alloc,'r') as file:
    alloc_info = file.readlines()
    for i in range(1,len(alloc_info)):
        id_info = id_info+[int(alloc_info[i].split()[0])]
        ix      = ix+[int(alloc_info[i].split()[8])-1]
        iy      = iy+[int(alloc_info[i].split()[9])-1]
id_info  = np.array(id_info)
ilon     = np.array(ix)
ilat     = np.array(iy)

def cal_situ_ref(id_info_new):
    ref_new = np.full(np.size(id_info_new),np.nan)
    for i in range(0,np.size(id_info_new)):
        where_ref = np.where(site_name==id_info_new[i])
        ref_new[i]= ref[where_ref]
    ref_new = np.where((ref_new>998.99)&(ref_new<999.5),np.nan,ref_new)
    return ref_new
site_ref = cal_situ_ref(id_info)


def read_nc(file,num):
    nf = nc.Dataset(file,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    print(varname)
    lat = np.array(nf.variables[varname[0]][:])
    lon = np.array(nf.variables[varname[1]][:])
    var = np.array(nf.variables[varname[num]][:])
    var = np.where(var<-1000,np.nan,var)
    var = np.where(var>10**8,np.nan,var)
    return lat,lon,var

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

#%% read data
map_path = '/work/a06/yingying/CaMa_v411/cmf_v411_pkg/map/tej_01min/'
def read_bin(map_path,filename):
    with open(map_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
        day_array = day_array.reshape(1320,1500)
        day_array = np.where(day_array<-999,np.nan,day_array)
    return day_array

def read_nextxy(filename):
    with open(map_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.int32)
        day_array = day_array.reshape(2,1320,1500)
        day_array = np.where(day_array<-999,np.nan,day_array)
    return day_array

#%% read ens range
typeA = "A" # online DA
typeC = "C" # simulation
typeF = "F" # offline DA
output_path  = '/work/a06/yingying/Code/plot/'+exp+dahour_str+'/'

if var_name == 'outflw':
    lat,lon,disC_max = read_ens_nc(output_path+var_dis+'/data'+typeC+'_max01_20.nc')
    lat,lon,disC_min = read_ens_nc(output_path+var_dis+'/data'+typeC+'_min01_20.nc')
    lat,lon,disC_mean = read_ens_nc(output_path+var_dis+'/data'+typeC+'_mean01_20.nc')
    #lat,lon,disA_max = read_ens_nc(output_path+var_dis+'/data'+typeA+'_max01_20.nc')
    #lat,lon,disA_min = read_ens_nc(output_path+var_dis+'/data'+typeA+'_min01_20.nc')
    lat,lon,disA_mean = read_ens_nc(output_path+var_dis+'/data'+typeA+'_mean01_20.nc')
else:
    lat,lon,disC_max = read_ens_nc(output_path+var_dep+'/data'+typeC+'_max01_20.nc')
    lat,lon,disC_min = read_ens_nc(output_path+var_dep+'/data'+typeC+'_min01_20.nc')
    lat,lon,disC_mean = read_ens_nc(output_path+var_dep+'/data'+typeC+'_mean01_20.nc')
    #lat,lon,disA_max = read_ens_nc(output_path+var_dep+'/data'+typeA+'_max01_20.nc')
    #lat,lon,disA_min = read_ens_nc(output_path+var_dep+'/data'+typeA+'_min01_20.nc')
    lat,lon,disA_mean = read_ens_nc(output_path+var_dep+'/data'+typeA+'_mean01_20.nc')

basin=read_bin(map_path,'basin.bin')
up_area = read_bin(map_path,'uparea.bin')
elvmean = read_bin(map_path,'elv_mean.bin')
elvmean = np.where(elvmean<-100,np.nan,elvmean)
nextlon = read_nextxy('nextxy.bin')[0,:,:]-1
nextlat = read_nextxy('nextxy.bin')[1,:,:]-1
nextlon = np.where(nextlon<0,np.nan,nextlon)
nextlat = np.where(nextlat<0,np.nan,nextlat)

# read offline DA result
lat,lon,dataF = read_nc(output_path+var_dep+'/dataF.nc',3)    

#%% read rainfall data
def read_prep(file):
    nf = nc.Dataset(file,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    print(varname)
    lat = np.array(nf.variables[varname[0]][:])
    lon = np.array(nf.variables[varname[1]][:])
    var = np.array(nf.variables[varname[3]][:])
    var = var*3600 # unit: mm/h 
    var = np.array(np.where(var<-900,np.nan,var))
    return lat,lon,var
lat,lon,rainf = read_prep('/work/a06/yingying/ils_data/tej_rain.nc')

def write_bin(var_each,varname,num):
    var_flat = var_each.flatten()
    var_bin  = var_flat.astype('int32').tobytes()
    # write .bin
    fname    = "basin"+varname+str(num)
    filebin  = map_path+fname+'.bin'
    with open(filebin,'wb') as file:
        file.write(var_bin)
    
obs_path = '/work/a06/yingying/MLIT/data/data_'+var_dir+'_hour2019/'
def read_obs(situ_id):
    filename = str(situ_id)+'.bin'
    if os.path.exists(filename):            
        with open(obs_path+filename, 'rb') as f:
            day_array = np.fromfile(f, dtype=np.float32)
            day_array = np.where(day_array<-900,np.nan,day_array)
        f.close()
        return day_array
    else:
        return -999

#%%  draw the da output
os.makedirs('/work/a06/yingying/plot/'+exp+dahour_str+'/',exist_ok=True)
os.makedirs('/work/a06/yingying/plot/'+exp+dahour_str+'/'+var_dir+'/',exist_ok=True)
def draw_lines(loc_lat,loc_lon,ylabel,loc_sim_mean,loc_sim_max,loc_sim_min,loc_da,loc_obs,locF,rain):
    fig,ax1 = plt.subplots(dpi = 300,figsize=(6,3))
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
    time    = np.arange(0,time_len,1)
    loc_lat = '{:04d}'.format(loc_lat)
    loc_lon = '{:04d}'.format(loc_lon)   
    # dis obs
    ax1.plot(time,loc_obs,color='blue',label='Obs',linewidth=1.5)
    # simulations
    ax1.plot(time,loc_sim_mean,'r-',linewidth=1.2,label='Sim')
    ax1.fill_between(time, loc_sim_max, loc_sim_min, color='red', alpha=0.2)
    # offline DA
    for day in range(0,int((day2-day1+1)*24/dahour)-1):
        time_each = np.arange(day*dahour,(day+1)*dahour,1)
        if var_name == "rivdph":
            ax1.scatter(time_each[-1],locF[day], color='black', marker='*', s=30)
        sday = day*dahour
        eday = (day+1)*dahour
        ax1.plot(time_each,loc_da[sday:eday],'g-',linewidth=1.2,alpha=0.7)       
        time_each_da =np.arange(sday,eday,1)
        ax1.scatter(time_each_da[-1],loc_da[eday-1],marker='.',color='green',alpha=0.9,s=30)
    # online DA
    day = day+1
    time_each = np.arange(day*dahour,(day+1)*dahour,1)
    ax1.plot(time_each,loc_da[day*dahour:(day+1)*dahour],'g-',linewidth=1.2,alpha=0.7,label='DA')
    if var_name == 'rivdph':
        ax1.scatter(time_each[-1],locF[-1], color='black', marker='*', s=30,label='Offline DA')

    ax1.set_xlabel('Hour')
    ax1.set_ylabel(ylabel)
    ax1.set_xlim(0,time_len)
    ax1.set_title('DA with Anomaly Water Depth '+str(loc_lon)+str(loc_lat))    
    ax1.legend(loc='upper left')

    dates = pd.date_range(start='2019-10-01 09:00',end='2019-10-19 08:00',freq='h')
    ax1.set_xticks(np.arange(15,time_len,24*3))
    ax1.set_xticklabels(dates[15::24*3],rotation=15)

    ax2 = ax1.twinx()
    ax2.invert_yaxis()
    ax2.bar(time,rain,color='gray',alpha=0.8)
    ax2.set_ylabel('Prep (mm/h)')
    ax2.set_ylim(np.nanmax(rain)+25,0)
    ax2.spines['right'].set_color('gray')
    ax2.yaxis.label.set_color('gray')
    ax2.tick_params(axis='y',colors='gray')
    # plt.show()
    plt.savefig('/work/a06/yingying/plot/'+exp+dahour_str+'/'+var_dir+'/'+str(loc_lon)+str(loc_lat)+'.jpg', format='jpg',dpi=300)
    plt.close()

def draw_obs(obs_lat,obs_lon):
    simmean = np.nanmean(disC_mean,axis=0)
    for loc in range(0,np.shape(obs_lat)[0]):
        loc_lat = obs_lat[loc]
        loc_lon = obs_lon[loc]
        loc_ref = site_ref[loc]
        loc_id  = id_info[loc]
        #obs
        obs_path = '/work/a06/yingying/MLIT/data/data_'+var_dir+'_hour2019/'
        filename = obs_path+str(loc_id)+'.bin'
        if os.path.exists(filename):            
            with open(filename, 'rb') as f:
                obs = np.fromfile(f, dtype=np.float32)
                obs = np.where(day_array<-900,np.nan,day_array)
            # simulation
            if var_name == 'outflw':
                # Japanese time is 9 hours later than UTC
                disC_sim_mean = disC_mean[0:time_len,loc_lat,loc_lon]
                disC_sim_max = disC_max[0:time_len,loc_lat,loc_lon]
                disC_sim_min = disC_min[0:time_len,loc_lat,loc_lon]
                disA_da_mean = disA_mean[0:time_len,loc_lat,loc_lon]
                dis_obs = obs[0:time_len]
                disF_loc= np.full((int((day2-day1+1)*24/dahour),1320,1500),np.nan)
            else:
                disC_sim_mean = disC_mean[0:time_len,loc_lat,loc_lon]-simmean[loc_lat,loc_lon]
                disC_sim_max = disC_max[0:time_len,loc_lat,loc_lon]-simmean[loc_lat,loc_lon]
                disC_sim_min = disC_min[0:time_len,loc_lat,loc_lon]-simmean[loc_lat,loc_lon]
                disA_da_mean = disA_mean[0:time_len]-simmean[loc_lat,loc_lon]
                dis_obs = obs[0:time_len,loc_lat,loc_lon]-elvmean[loc_lat,loc_lon]
                #offline DA
                disF_loc= dataF[:,loc_lat,loc_lon]
            if np.all(disA_mean[50:100,loc_lat,loc_lon]<10**(-8))==True:
                disA_da_mean = np.full(time_len,np.nan)
            rain    = rainf[0:time_len,loc_lat,loc_lon]
            draw_lines(loc_lat,loc_lon,var_ylabel,disC_sim_mean,disC_sim_max,disC_sim_min,disA_da_mean,dis_obs,disF_loc,rain)
        else:
            print(var_name+': no observation')

#draw_lines(loc_lat,loc_lon,var_ylabel,disC_sim_mean,disC_sim_mean,disC_sim_mean,disA_da_mean,dis_obs,disF_loc,rain)

# draw_obs([564-1],[1006-1])
# draw_obs([421-1],[1095-1])
draw_obs(ilat,ilon)

#%% NRMSE
def draw_each_figure(obs_lat,obs_lon):
    situ_rmse  = list()
    situ_nrmse  = list()
    lat_plot   = list()
    lon_plot   = list()
    situ_imp   = list()
    simmean = np.nanmean(disC_mean,axis=0)
    
    for loc in range(0,np.shape(obs_lat)[0]):
        loc_lat = obs_lat[loc]
        loc_lon = obs_lon[loc]
    # read data
        if var_name == 'outflw':
            loc_sim = disC_mean[0:time_len,loc_lat,loc_lon]
            loc_obs = obs[0:time_len,loc_lat,loc_lon]
            loc_da  = disA_mean[0:time_len,loc_lat,loc_lon]
        else:
            loc_sim = disC_mean[0:time_len,loc_lat,loc_lon]-simmean[loc_lat,loc_lon]
            loc_obs = obs[0:time_len,loc_lat,loc_lon]-elvmean[loc_lat,loc_lon]
            loc_da  = disA_mean[0:time_len,loc_lat,loc_lon]-simmean[loc_lat,loc_lon]
        if np.all(disA_mean[50:100,loc_lat,loc_lon]<10**(-8))==True:
            loc_da = np.full(time_len,np.nan)
        # remove no obs stations   
        if np.all(np.isnan(loc_obs))==True | np.all(loc_obs==0) | np.all(np.isnan(loc_sim))==True | np.all(np.isnan(loc_da))==True:
            continue
        where_nan = np.where(np.isnan(loc_obs))[0]
        if len(where_nan>0):
            loc_obs   = np.delete(loc_obs,where_nan)
            loc_sim   = np.delete(loc_sim,where_nan)
            loc_da    = np.delete(loc_da,where_nan)
            # print(where_nan)
        where_nan = np.where(np.isnan(loc_da))[0]
        if len(where_nan>0):
            loc_obs   = np.delete(loc_obs,where_nan)
            loc_sim   = np.delete(loc_sim,where_nan)
            loc_da    = np.delete(loc_da,where_nan)
        where_nan = np.where(np.isnan(loc_sim))[0]
        if len(where_nan>0):
            loc_obs   = np.delete(loc_obs,where_nan)
            loc_sim   = np.delete(loc_sim,where_nan)
            loc_da    = np.delete(loc_da,where_nan)
        if (np.shape(loc_obs)[0]<5) | (np.shape(loc_sim)[0]<5) | (np.shape(loc_da)[0]<5):
            rmse = np.nan
            nrmse = np.nan
            imp = np.nan
        else:
            mse_sim  = mean_squared_error(loc_obs, loc_sim)
            mse_da   = mean_squared_error(loc_obs, loc_da)
            rmse     = np.sqrt(mse_sim)-np.sqrt(mse_da)
            var_max  = np.nanmax(loc_obs)
            var_min  = np.nanmin(loc_obs)
            var_range = var_max-var_min
            nrmse    = rmse/var_range
            if nrmse<-0.2:
                print(nrmse,rmse,loc_lat,loc_lon)
            imp = rmse/np.sqrt(mse_sim)
        # store nrmse, lat, loc
        situ_rmse.append(rmse)
        situ_nrmse.append(nrmse)
        situ_imp.append(imp)
        lat_plot.append(loc_lat)
        lon_plot.append(loc_lon)
    situ_rmse = np.array(situ_rmse)
    situ_nrmse = np.array(situ_nrmse)
    lat_plot = np.array(lat_plot)
    lon_plot = np.array(lon_plot)
    situ_imp = np.array(situ_imp)
    return situ_rmse,situ_nrmse,lat_plot,lon_plot,situ_imp
situ_rmse,situ_nrmse,lat_plot,lon_plot,situ_imp = draw_each_figure(ilat,ilon)

nrmse_pos = np.shape(np.where(situ_nrmse>0))[1]/(np.shape(np.where(situ_nrmse>0))[1]+np.shape(np.where(situ_nrmse<0))[1])
nrmse_neg = 1-nrmse_pos
print('better',nrmse_pos,'worse',nrmse_neg)
print('average performance',np.nanmean(situ_nrmse),np.nanmean(situ_imp))
#%% draw the map
def draw_rmse(title,situ_nrmse,obs_lon,obs_lat,vmin,vmax):
    fig = plt.figure(dpi = 600,figsize=(7.5,4.5))
    ax1 = fig.add_axes([0.1,0.1,0.8,0.8],projection=ccrs.PlateCarree())
    # ax1.outline_patch.set_linewidth(0.35)
    ax1.set_extent([128,149,30,46], ccrs.PlateCarree())
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=1, color='gray', alpha=0, linestyle='-.')
    ax1.coastlines(alpha=1.,linestyle='-',color = "#3F3A3A",lw=0.8,resolution='10m')
    
    # mask the other place
    region_lon1 = [128.06, 128.06,138, 138 ]
    region_lat1 = [41 , 45.95 , 45.95, 41  ]
    region_lon2 = [128.06,128.06,131.4,131.4 ]
    region_lat2 = [34.5,41,41,34.5]
    ax1.fill(region_lon1, region_lat1, color='white', transform=ccrs.PlateCarree(), zorder=10)
    ax1.fill(region_lon2, region_lat2, color='white', transform=ccrs.PlateCarree(), zorder=10)
    
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([130,135,140,145])
    gl.ylocator = mticker.FixedLocator([30,35,40,45])
    gl.xlabel_style = {'size': 8,'color':"#3F3A3A",'weight':'normal'}
    gl.ylabel_style = {'size': 8,'color':"#3F3A3A",'weight':'normal'}
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)
    norm = mcolors.TwoSlopeNorm(vmin=min(situ_rmse), vcenter=0, vmax=-min(situ_rmse))
    con_plot2 = ax1.scatter(lon[obs_lon],lat[obs_lat],c=situ_nrmse,s=1.5,zorder=2,cmap='coolwarm',transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,alpha=0.8)

    def color_bar(l,b,w,h):
      rect = [l,b,w,h]
      cbar_ax = fig.add_axes(rect)
      return cbar_ax
    [l1,b1,w1,h1] = [0.63,0.12,0.01,0.25]
    [l2,b2,w2,h2] = [0.73,0.12,0.01,0.25]
    cbar_ax2 = color_bar(l2,b2,w2,h2)

    # cb1.ax.set_title('NRMSE \n $(m^{3}/s)$',fontdict=font_label)
    cb2 = plt.colorbar(con_plot2, cax=cbar_ax2,orientation="vertical",shrink=0.5)
    ax1.set_title(title,loc='center')
    # plt.show()
    plt.savefig('/work/a06/yingying/plot/'+exp+dahour_str+'/'+title+'.jpg', format='jpg',dpi=600)
    # plt.close()
if var_name == 'rivdph':
    draw_rmse(var_name+'_RMSE',situ_rmse,lon_plot,lat_plot,-1,1)   
else:
    draw_rmse(var_name+'_RMSE',situ_rmse,lon_plot,lat_plot,-max(situ_rmse)*0.6,max(situ_rmse)*0.6)   
draw_rmse(var_name+'_nRMSE',situ_nrmse,lon_plot,lat_plot,-0.4,0.4)   
draw_rmse(var_name+'_Improvement',situ_imp,lon_plot,lat_plot,0,1)   
