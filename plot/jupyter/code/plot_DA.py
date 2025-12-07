import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case21 as pm
pm_path = pm.__file__
pm_name = pm.__name__

import numpy as np
import netCDF4 as nc
import glob
import os
import pandas as pd
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,BoundaryNorm
import matplotlib.cm as cm
from matplotlib import colors
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import folium
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.io import shapereader

def cal_assimilate_grid():
    with open(pm.CaMa_dir()+'/map/tej_01min/nextxy.bin', 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.int32)
        data_re  = np.reshape(data_bin,(2,1320,1500))
    flow_ix = data_re[0,:,:]-1  # next ix (ilon)
    flow_iy = data_re[1,:,:]-1  # next iy (ilat)
    flow_ix = np.where((flow_ix<0),np.nan,flow_ix)
    flow_iy = np.where((flow_iy<0),np.nan,flow_iy)
    patchdir = '/data50/yingying/HydroDA/obs_patch/obs_patch06/obs_patch/'
    txt_files= glob.glob(patchdir+'*.txt')
    patchname = [os.path.splitext(f)[0][-8:] for f in txt_files]
    data_weight = np.load('/data50/yingying/HydroDA/src/plot/weight.npy')
    with open(pm.CaMa_dir()+'/map/tej_01min/uparea.bin', 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(1320,1500))
    masked_data = np.full_like(data_re, np.nan)  
    valid = (~np.isnan(flow_ix)) & (~np.isnan(flow_iy))
    ix = flow_ix[valid].astype(int)
    iy = flow_iy[valid].astype(int)
    masked_data[iy, ix] = data_re[iy, ix]
    masked_data = np.where(masked_data<10**8,np.nan,masked_data)  # 100m2
    print(np.shape(np.where(~np.isnan(masked_data))))
    river_grid_num = np.shape(np.where(~np.isnan(masked_data)))[1]
    patch_num = np.shape(data_weight)[0]
    rate = patch_num/river_grid_num*100
    rateg = 1500/river_grid_num*100
    print('assimilated rate:','{:.2f}'.format(rate),'gauge rate:','{:.2f}'.format(rateg))    
    ## [lat, lon, weight]
    # data_weight = np.full((patch_num,3),np.nan)
    # for ind in range(len(txt_files)):
    #     data  = np.array(np.loadtxt(txt_files[ind]))
    #     if data.ndim == 2:
    #         data_weight[ind,2] = np.nanmax(data[:,2])
    #     else:
    #         data_weight[ind,2] = data[2]
    #     data_weight[ind,1] = int(patchname[ind][:4])-1
    #     data_weight[ind,0] = int(patchname[ind][4:])-1
    return data_weight,masked_data#,lon_flow,lat_flow

#%% draw the select stations
def draw_weight(data_weight,lon_obs,lat_obs,lon,lat,data_uparea,inputdir):
    fig = plt.figure(dpi = 600,figsize=(7.5,4.5))
    ax1 = fig.add_axes([0.1,0.1,0.8,0.8],projection=ccrs.PlateCarree())
    ax1.set_axis_off()
    # ax1.outline_patch.set_linewidth(0.35)
    ax1.set_extent([128,149,30,46], ccrs.PlateCarree())
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=1, color='gray', alpha=0, linestyle='-.')
    ax1.coastlines(alpha=1.,linestyle='-',color = "#3F3A3A",lw=0.6,resolution='10m')
    
    # mask the other place
    region_lon1 = [128.06, 128.06,138, 138 ]
    region_lat1 = [41 , 45.95 , 45.95, 41  ]
    region_lon2 = [128.06,128.06,131.4,131.4 ]
    region_lat2 = [34.5,41,41,34.5]
    ax1.fill(region_lon1, region_lat1, color='white', transform=ccrs.PlateCarree(), zorder=10)
    ax1.fill(region_lon2, region_lat2, color='white', transform=ccrs.PlateCarree(), zorder=10)
    
    gl.top_labels   = False
    gl.right_labels = False
    gl.bottom_labels   = False
    gl.left_labels = False
    lon_weight_ind = np.array(data_weight[:,1],np.int32)
    lat_weight_ind = np.array(data_weight[:,0],np.int32)
    cmap_blue = ListedColormap(['gray']) 
    d_lon = lon[1] - lon[0]
    d_lat = lat[1] - lat[0]
    lon_edge = np.append(lon - d_lon / 2, lon[-1] + d_lon / 2)
    lat_edge = np.append(lat - d_lat / 2, lat[-1] + d_lat / 2)
    ax1.pcolormesh(lon_edge,lat_edge,data_uparea,cmap=cmap_blue,zorder=1, transform=ccrs.PlateCarree(),alpha=1.)  
    size=4
    # ax1.pcolormesh(lon2d_low_map,lat2d_low_map,up_low,cmap=cmap_blue,zorder=1, transform=ccrs.PlateCarree(),alpha=1.)
    con_plot1 = ax1.scatter(lon_obs,lat_obs,s=5,zorder=2,edgecolor='#3F3A3A',facecolors='none',lw=0.35,transform=ccrs.PlateCarree())
    con_plot2 = ax1.scatter(lon[lon_weight_ind],lat[lat_weight_ind],c=data_weight[:,2],cmap='rainbow',s=0.003,zorder=3,lw=1.,transform=ccrs.PlateCarree())
    def color_bar(l,b,w,h):
      rect = [l,b,w,h]
      cbar_ax = fig.add_axes(rect)
      return cbar_ax
    [l1,b1,w1,h1] = [0.75,0.12,0.01,0.25]
    [l2,b2,w2,h2] = [0.75,0.32,0.01,0.25]
    cbar_ax2 = color_bar(l2,b2,w2,h2)
    cb_title = 'localization parameter'
    cb2 = plt.colorbar(con_plot2, cax=cbar_ax2,orientation="vertical",shrink=0.5)
    cb2.ax.set_title(cb_title,fontsize=10)
    
    plt.savefig(inputdir+'/1.weight.jpg', format='jpg',dpi=600)
    plt.show()
    plt.close()

# plot relationship between observation numbers and RMSE, KGE
def plot_obs_relationship(ax,rmse_wlv):
    assim_uni = np.unique(rmse_wlv[:,7])
    assim_uni = np.array(assim_uni,dtype=np.int32)
    rmse_assim_da = np.full(len(assim_uni),np.nan)
    rmse_assim_sim = np.full(len(assim_uni),np.nan)
    med = list()
    q25 = list()
    q75 = list()
    q10 = list()
    q90 = list()
    for ind in range(0,len(assim_uni)):
        # rmse_sim-rmse_da (>0: good; <0: bad)
        rmse_assim = rmse_wlv[:,2]-rmse_wlv[:,1]
        med = med + [np.nanpercentile(rmse_assim[rmse_wlv[:,7]==assim_uni[ind]], 50)]
        q25 = q25 + [np.nanpercentile(rmse_assim[rmse_wlv[:,7]==assim_uni[ind]], 25)]
        q75 = q75 + [np.nanpercentile(rmse_assim[rmse_wlv[:,7]==assim_uni[ind]], 75)]
        q10 = q10 + [np.nanpercentile(rmse_assim[rmse_wlv[:,7]==assim_uni[ind]], 10)]
        q90 = q90 + [np.nanpercentile(rmse_assim[rmse_wlv[:,7]==assim_uni[ind]], 90)]
    err_lower = np.array(med) - np.array(q10)
    err_upper = np.array(q90) - np.array(med)
    x = np.arange(len(assim_uni))
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1.5)
    bar_width = 0.4
    errors = np.array([err_lower,err_upper])
    ax.errorbar(x, med, yerr=errors, fmt='o', capsize=5, color='black')
    ax.set_ylabel(r'RMSE$_{sim}$-RMSE$_{DA}$ (m)',fontsize=14)
    ax.set_xlabel('Observation numbers for assimilation', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(assim_uni)


# draw single stations for poster
def draw_single_station(i,lat_all_flood, lon_all_flood,obs_station_all, da_station_all, sim_station_all,syyyy,tstart,time_len,simmean,pm_name,rainf,lon,lat,pm_path,expname,dahour_str,exp_plotdir):
    time    = np.arange(0,time_len,1)
    time_plot= time[tstart:]
    fig,axes = plt.subplots(dpi = 600,figsize=(12,6))
    obs_plot = obs_station_all[tstart:,i]
    loc_lat  = lat_all_flood[i]-1
    loc_lon  = lon_all_flood[i]-1
    mean_station = simmean[int(loc_lat),int(loc_lon)]
    sim_plot = sim_station_all[tstart:,:,i]-mean_station
    da_plot  = da_station_all[tstart:,:,i]-mean_station 
    # prep
    ax2 = axes.twinx()
    ax2.invert_yaxis()
    ax2.bar(time[tstart:],rainf[tstart:,int(loc_lat),int(loc_lon)],color='gray',alpha=0.8)
    ax2.set_ylim(np.nanmax(rainf[tstart:,int(loc_lat),int(loc_lon)])+5,0)
    ax2.spines['right'].set_color('gray')
    ax2.yaxis.label.set_color('gray')
    ax2.tick_params(axis='y',colors='gray')
    # obs
    axes.plot(time_plot,obs_plot,color='green',linewidth=2.,linestyle='--')
    # sim
    axes.plot(time_plot,sim_plot[:,-1],'b-',linewidth=3.5)
    # DA
    for size in range(np.shape(da_plot)[1]):
        axes.plot(time_plot,da_plot,'r-',linewidth=1.25)
    axes.set_xlim(tstart,time_len)
    if 'case' in pm_name: # parmas_case11
        if syyyy == '2022':
            dates = pd.date_range(start='2022-08-03 00:00',end='2022-08-06 08:00',freq='24h')
            axes.set_xticks(np.arange(tstart,time_len,24))
        else: # syyyy='2024'
            dates = pd.date_range(start='2024-07-24 08:00',end='2024-07-25 07:00',freq='8h')
            axes.set_xticks(np.arange(tstart,time_len,8))
    else:  # syyyy ='2019'
        dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='48h')
        axes.set_xticks(np.arange(tstart,time_len,48))
    axes.set_xticklabels(dates,fontsize=15)
    axes.tick_params(axis='y',labelsize=15)
    ax2.tick_params(axis='y',labelsize=15)
    plt.savefig(f'./agu/{syyyy}.jpg',dpi=600)
    plt.show()
    plt.close()
    
