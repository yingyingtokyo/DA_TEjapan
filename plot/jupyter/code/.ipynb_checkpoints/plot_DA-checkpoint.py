import numpy as np
import netCDF4 as nc
import glob
import os

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

import importlib
import params_ils01 as pm    # out_msm_round*  # params_all_case3
pm_path = pm.__file__
pm_name = pm.__name__

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
def draw_weight(data_weight,lon_obs,lat_obs,lon,lat,inputdir):
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