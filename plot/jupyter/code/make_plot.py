import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case21 as pm
pm_path = pm.__file__
pm_name = pm.__name__

import numpy as np
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
import matplotlib.ticker as mticker
import pandas as pd

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.io import shapereader

import importlib
import basic_command as basic
importlib.reload(basic)
from basic_command import find_id, read_bin2d, read_bin3d, read_bin_int, read_bin_station, read_nc,cal_rmse,cal_kge,cal_nse,read_bin_mean

#%% draw the select stations
def draw_select_station(lon_valid,lat_valid,assim_num,lon_select,lat_select,lon2d_low_map,lat2d_low_map,up_low,inputdir,pm_path,title):
    fig = plt.figure(dpi = 600,figsize=(7.5,4.5))
    ax1 = fig.add_axes([0.1,0.1,0.8,0.8],projection=ccrs.PlateCarree())
    ax1.set_axis_off()
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
    gl.bottom_labels   = False
    gl.left_labels = False
    # con_plot1 = ax1.scatter(lon[ix],lat[iy],s=1,zorder=2,color='black',transform=ccrs.PlateCarree(),alpha=0.95)
    con_plot1 = ax1.scatter(lon_select,lat_select,s=1,zorder=2,color='gray',transform=ccrs.PlateCarree(),alpha=0.95)
    # if pm.val_mode()==6:
    #     assim_num = 2
    con_plot2 = ax1.scatter(lon_valid,lat_valid,s=assim_num*2,zorder=3,edgecolor='red',facecolors='none',lw=1.,transform=ccrs.PlateCarree())
    def color_bar(l,b,w,h):
      rect = [l,b,w,h]
      return cbar_ax
    [l1,b1,w1,h1] = [0.63,0.12,0.01,0.25]
    [l2,b2,w2,h2] = [0.63,0.32,0.01,0.25]
    # cbar_ax2 = color_bar(l2,b2,w2,h2)
    cmap_blue = ListedColormap(['skyblue']) 
    ax1.pcolormesh(lon2d_low_map,lat2d_low_map,up_low,cmap=cmap_blue,zorder=1, transform=ccrs.PlateCarree(),alpha=1.)  

    # cb1.ax.set_title('NRMSE \n $(m^{3}/s)$',fontdict=font_label)
    # cb2 = plt.colorbar(con_plot2, cax=cbar_ax2,orientation="vertical",shrink=0.5)
    # ax1.set_title(title,loc='center')
    if 'round' in pm_path:
        plt.savefig(inputdir+title+'_round'+pm_path[-4]+'.jpg', format='jpg',dpi=600)
    elif 'river' in pm_path:
        plt.savefig(inputdir+title+'_river'+pm_path[-4]+'.jpg', format='jpg',dpi=600)
    elif 'case' in pm_path:
        plt.savefig(inputdir+title+'_case'+pm_path[-5]+pm_path[-4]+'.jpg', format='jpg',dpi=600)
    else:
        plt.savefig(inputdir+title+'.jpg', format='jpg',dpi=600)
    plt.show()
    plt.close()


def draw_lines(varname,lat_valid_plot, lon_valid_plot, assim_num_plot, obs_station,da_station,sim_station,savename,new_row,alabel,syyyy,tstart,time_len,simmean,pm_name,rainf,lon,lat,pm_path,expname,dahour_str,exp_plotdir):
    if 'case' in pm_name:
        if syyyy == '2022':
            rows_old, cols = 11, 3
            fig,axes = plt.subplots(rows_old, cols, dpi = 600,figsize=(15,28))
        else:
            rows_old, cols = 6, 3
            fig,axes = plt.subplots(rows_old, cols, dpi = 600,figsize=(15,24))
        obs_range=0.01
    else:
        rows_old, cols = 10, 3
        if "WSE" in varname:
            fig,axes = plt.subplots(rows_old, cols, dpi = 600,figsize=(15,35))
            obs_range=0.01
        else:
            fig,axes = plt.subplots(rows_old, cols, dpi = 600,figsize=(15,28))
            obs_range=10
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
    handles = [
    Line2D([0], [0], color='green', lw=2, label='In-situ observation',linestyle='--'),
    Line2D([0], [0], color='blue', lw=2, label='Open-loop result without rainfall pertubation'),
    Line2D([0], [0], color='blue', lw=0.4, label='Open-loop ensembles'),
    Line2D([0], [0], color='red', lw=2, label='DA ensembles')]

    fig.legend(handles=handles,fontsize=12,
        loc='upper center',
        ncol=4,
        bbox_to_anchor=(0.5, 0.925), 
        frameon=False)
    plt.subplots_adjust(bottom=0.15) 
    time = np.arange(0,time_len,1)
    valid_indices = [i for i in range(np.shape(obs_station)[1])]
    improve_rate = np.zeros(len(valid_indices))
    num_valid = len(valid_indices)
    rows = (num_valid + cols - 1) // cols
    # plot stations
    plot_idx = 0
    lon_point = []
    lat_point = []
    lon_point_ind = []
    lat_point_ind = []
    obs_point_ind = []
    sim_point_ind = []
    da_point_ind = []
    num_point = []
    kge_point = []
    nse_point = []
    # [plot_idx, rmse_da, rmse_sim, kge_da, kge_sim, lat, lon, assim_num]
    rmse_all = []
    for ind, i in enumerate(valid_indices):
        time_plot= time[tstart:]
        obs_plot = obs_station[tstart:,i]
        loc_lat  = lat_valid_plot[i]-1
        loc_lon  = lon_valid_plot[i]-1
        mean_station = simmean[int(loc_lat),int(loc_lon)]
        sim_plot = sim_station[tstart:,:,i]-mean_station
        da_plot  = da_station[tstart:,:,i]-mean_station
        # RMSE
        rmse_da = cal_rmse(obs_plot,np.nanmean(da_plot,axis=1),varname,time_len)
        kge_da = cal_kge(obs_plot,np.nanmean(da_plot,axis=1),varname)
        nse_da = cal_nse(obs_plot,np.nanmean(da_plot,axis=1),varname)
        rmse_da_str = '{:.2f}'.format(rmse_da)
        rmse_sim= cal_rmse(obs_plot,sim_plot[:,-1],varname,time_len)  # compare with open loop
        kge_sim= cal_kge(obs_plot,sim_plot[:,-1],varname)  # compare with open loop
        nse_sim= cal_nse(obs_plot,sim_plot[:,-1],varname)  # compare with open loop
        rmse_sim_str= '{:.2f}'.format(rmse_sim)
        nse_point = nse_point + [[nse_da,nse_sim]]
        kge_point = kge_point + [[kge_da,kge_sim]]
        if (np.isnan(rmse_da)) | (np.isnan(rmse_sim)):
            improve_rate[ind] = np.nan
            continue
        if (np.shape(np.unique(da_plot))[0]<2):  # exclude the fail DA result
            improve_rate[ind] = np.nan
            continue
        if 'case1' in pm_name:
            if (np.abs(np.nanmax(obs_plot))<0.5):  # exclude the gauges at small rivers
                improve_rate[ind] = np.nan
                continue 
        else:
            if (np.abs(np.nanmax(obs_plot))<0.1):  # exclude the gauges at small rivers
                improve_rate[ind] = np.nan
                continue            
        # Rainfall
        row, col = divmod(plot_idx,cols)
        ax2 = axes[row,col].twinx()
        ax2.invert_yaxis()
        ax2.bar(time[tstart:],rainf[tstart:,int(loc_lat),int(loc_lon)],color='gray',alpha=0.8)
        ax2.set_ylim(np.nanmax(rainf[tstart:,int(loc_lat),int(loc_lon)])+5,0)
        ax2.spines['right'].set_color('gray')
        ax2.yaxis.label.set_color('gray')
        ax2.tick_params(axis='y',colors='gray') 
        # simulations
        sim_point_ind = sim_point_ind + [sim_plot]
        axes[row,col].plot(time_plot,sim_plot,'b-',linewidth=0.2,alpha=0.8)
        # online DA
        da_point_ind  = da_point_ind + [da_plot]
        axes[row,col].plot(time_plot,da_plot,'r-',linewidth=1.25)
        axes[row,col].set_xticks(np.arange(tstart+15,time_len,24))
        axes[row,col].set_xlim(tstart,time_len)
        # open loop
        axes[row,col].plot(time_plot,sim_plot[:,-1],'b-',linewidth=1.5)
        # axes[row,col].plot(time_plot,np.nanmean(sim_plot,axis=1),'b-',linewidth=1.5)
        # obs
        obs_point_ind = obs_point_ind + [obs_plot]
        axes[row,col].plot(time_plot,obs_plot,color='green',linewidth=2.5,linestyle='--')
        # lat_str = '{:.2f}'.format(lat[int(lat_valid_plot[i])-1])
        # lon_str = '{:.2f}'.format(lon[int(lon_valid_plot[i])-1])
        lat_str = '{:.2f}'.format(int(lat_valid_plot[i])-1)
        lon_str = '{:.2f}'.format(int(lon_valid_plot[i])-1)
        if rmse_da < rmse_sim:
            improve_rate[ind] = 1
        rmse_all = rmse_all+ [plot_idx,rmse_da,rmse_sim,kge_da,kge_sim,int(lat_valid_plot[i])-1,int(lon_valid_plot[i])-1,int(assim_num_plot[i])]
        # axes[row,col].set_title(alabel[plot_idx],loc='left',fontsize=12)
        axes[row,col].text(0.02,0.95,'('+lat_str+','+lon_str+')\n' + 'obs number: '+str(int(assim_num_plot[i])),transform=axes[row,col].transAxes,fontsize=12, c='black',fontweight='bold', va='top', ha='left')
        axes[row,col].text(0.9,0.95,'RMSE: '+rmse_da_str,transform=axes[row,col].transAxes,fontsize=12, c='red',fontweight='bold', va='top', ha='right')
        axes[row,col].text(0.9,0.85,'RMSE: '+rmse_sim_str,transform=axes[row,col].transAxes,fontsize=12, c='blue',fontweight='bold', va='top', ha='right')
        if col == 0:
            axes[row,col].set_ylabel(varname)
        if col == cols-1:
            ax2.set_ylabel('Rain (mm/h)')
        if 'WSE' in varname:
            if 'case' in pm_name: # parmas_case11
                if syyyy == '2022':
                    dates = pd.date_range(start='2022-08-03 00:00',end='2022-08-06 08:00',freq='24h')
                else: # syyyy='2024'
                    dates = pd.date_range(start='2024-07-24 08:00',end='2024-07-25 07:00',freq='24h')
                axes[row,col].set_xticklabels(dates,rotation=25)
            else:  # syyyy ='2019'
                if row == 9: # parmas01
                    dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='24h')
                    axes[row,col].set_xticklabels(dates,rotation=25)
                else:
                    axes[row,col].set_xticklabels([],rotation=15)
        else:
                if row == 5:
                    dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='24h')
                    axes[row,col].set_xticklabels(dates,rotation=25)
                else:
                    axes[row,col].set_xticklabels([],rotation=25) 
        plot_idx = plot_idx + 1
        lon_point = lon_point + [lon[int(lon_valid_plot[i])-1]]
        lat_point = lat_point + [lat[int(lat_valid_plot[i])-1]]
        lon_point_ind  = lon_point_ind + [int(lon_valid_plot[i])-1]
        lat_point_ind  = lat_point_ind + [int(lat_valid_plot[i])-1]
        num_point = num_point + [assim_num_plot[i]]
    for ind in range(plot_idx, rows_old * cols):
        row, col = divmod(ind, cols)
        fig.delaxes(axes[row, col])  # delete the empty 
    improve = np.nansum(improve_rate)/np.shape(np.where(~np.isnan(improve_rate)))[1]*100
    print(np.nansum(improve_rate),np.shape(np.where(~np.isnan(improve_rate)))[1])
    print('How many stations improve: ','{:.2f}'.format(improve),'%')
    new_row = new_row + ['{:.2f}'.format(improve)]
    if 'all' in pm_path:
        # assimilated
        if 'case' in pm_path:
            plt.savefig(exp_plotdir+'2.cross_allcase'+savename+pm_path[-4]+'.jpg', format='jpg',dpi=600)
        else:
            plt.savefig(exp_plotdir+'2.cross_all'+savename+'.jpg', format='jpg',dpi=600)
    else:
        # validated
        if 'round' in pm_path:
            np.save(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse_round'+savename+pm_name[-1]+'.npy',np.array(rmse_all))
            plt.savefig(exp_plotdir+'2.cross_validation'+savename+pm_path[-4]+'.jpg', format='jpg',dpi=600)
        elif 'case' in pm_path:
            np.save(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse_case'+savename+pm_name[-2]+pm_name[-1]+'.npy',np.array(rmse_all))
            plt.savefig(exp_plotdir+'2.cross_valcase'+pm_name[-1]+savename+'.jpg', format='jpg',dpi=600)
        else:
            np.save(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse'+savename+'.npy',np.array(rmse_all))
            plt.savefig(exp_plotdir+'2.cross_validation'+savename+'.jpg', format='jpg',dpi=600)
    plt.show()
    plt.close()
    return np.array(lon_point), np.array(lat_point), np.array(lon_point_ind), np.array(lat_point_ind), np.array(obs_point_ind), np.array(sim_point_ind), np.array(da_point_ind), np.array(num_point), new_row


# plot station 24 hour forecast wlv validation
def plot_gradient_line(ax, x, y, cmap, linewidth=1.0, alpha=0.8, vmin=None, vmax=None,add_colorbar=False,cbar_label='Leading time',fig=None):
    x = np.asarray(x)
    y = np.asarray(y)
    if len(x) != len(y):
        raise ValueError("x and y must be the same length")
    if len(x) < 2:
        return  # Don't do anything
    if isinstance(cmap, str):
        cmap = plt.colormaps.get_cmap(cmap)
        colors = cmap(np.linspace(0.15, 0.85, len(x)-1))
    if vmin is None:
        vmin = 0
    if vmax is None:
        vmax = len(x) - 2
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    # check the results with different lead times
    # points = np.array([x[18:24], y[18:24]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=norm,
                        linewidth=linewidth, alpha=alpha)
    lc.set_array(np.arange(len(x)-1)) 
    ax.add_collection(lc)
    if add_colorbar and fig is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # no need for the real value
        ticks = np.linspace(0,len(x)-1,len(x)).astype(int)+1
        cbar = fig.colorbar(sm, ax=ax, pad=0.005)
        cbar.set_label(cbar_label,fontsize=15) 
        cbar.ax.tick_params(labelsize=15)  
        cbar.ax.yaxis.set_ticks_position('right') 
        cbar.set_ticks(ticks-1)
        cbar.set_ticklabels(ticks)
        cbar.ax.set_position([0.92, 0.4, 0.02, 0.7]) 
    return lc


 
def draw_fore_lines(varname,wlv_da,savename,station_lat,station_lon,station_idx,obs_station_all,lat_all_flood,lon_all_flood,simmean,sim_station_all,da_station_all,time_len,rainf,exp_plotdir,tstart,tstart_da,dahour,lat,lon):
    rows_old, cols = 14, 3
    if "WSE" in varname:
        fig,axes = plt.subplots(rows_old, cols, dpi = 300,figsize=(15,35))
        obs_range=0.01
    else:
        fig,axes = plt.subplots(rows_old, cols, dpi = 300,figsize=(15,28))
        obs_range=10
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
    handles = [
    Line2D([0], [0], color='green', lw=2, label='In-situ observation',linestyle='--'),
    Line2D([0], [0], color='blue', lw=2, label='Open-loop result without rainfall pertubation'),
    Line2D([0], [0], color='red', lw=0.5, label='DA forecasts')]

    fig.legend(handles=handles,fontsize=12,
        loc='upper center',
        ncol=3,
        bbox_to_anchor=(0.5, 0.925), 
        frameon=False)
    plt.subplots_adjust(bottom=0.15) 
    time    = np.arange(0,time_len,1)
    valid_indices = [i for i in range(len(lat_all_flood))]
    improve_rate = np.zeros(len(valid_indices))
    num_valid = len(valid_indices)
    rows = (num_valid + cols - 1) // cols
    # plot stations
    plot_idx = 0
    for ind, i in enumerate(valid_indices):
        time_plot= time[tstart:]
        obs_plot = obs_station_all[tstart:,i]
        loc_lat  = lat_all_flood[i]-1
        loc_lon  = lon_all_flood[i]-1
        mean_station = simmean[int(loc_lat),int(loc_lon)]
        sim_plot = sim_station_all[tstart:,:,i]-mean_station
        da_plot  = da_station_all[tstart:,:,i]-mean_station
        # RMSE
        rmse_da = cal_rmse(obs_plot,np.nanmean(da_plot,axis=1),varname,time_len)
        kge_da = cal_kge(obs_plot,np.nanmean(da_plot,axis=1),varname)
        nse_da = cal_nse(obs_plot,np.nanmean(da_plot,axis=1),varname)
        rmse_da_str = '{:.2f}'.format(rmse_da)
        rmse_sim= cal_rmse(obs_plot,sim_plot[:,-1],varname,time_len)  # compare with open loop
        if (np.isnan(rmse_da)) | (np.isnan(rmse_sim)):
            improve_rate[ind] = np.nan
            continue
        if (np.shape(np.unique(da_plot))[0]<2):  # exclude the fail DA result
            improve_rate[ind] = np.nan
            continue
        if (np.abs(np.nanmax(obs_plot))<0.1):  # exclude the gauges at small rivers
            improve_rate[ind] = np.nan
            continue
        # Rainfall
        row, col = divmod(plot_idx,cols)
        ax2 = axes[row,col].twinx()
        ax2.invert_yaxis()
        ax2.bar(time[tstart:],rainf[tstart:,int(loc_lat),int(loc_lon)],color='gray',alpha=0.8)
        ax2.set_ylim(np.nanmax(rainf[tstart:,int(loc_lat),int(loc_lon)])+5,0)
        ax2.spines['right'].set_color('gray')
        ax2.yaxis.label.set_color('gray')
        ax2.tick_params(axis='y',colors='gray') 
        # obs
        axes[row,col].plot(time_plot,obs_plot,color='green',linewidth=2.,linestyle='--')
        axes[row,col].set_xticks(np.arange(tstart,time_len,6))
        axes[row,col].set_xlim(tstart,time_len)
        axes[row,col].set_ylim(0,np.nanmax([np.nanmax(obs_plot),np.nanmax(np.nanmean(wlv_da[:,:,plot_idx,:],axis=0)),np.nanmax(sim_plot[:,-1])])+0.2)
        lat_str = '{:.2f}'.format(lat[int(loc_lat)])
        lon_str = '{:.2f}'.format(lon[int(loc_lon)])
        # use redish color to plot lines
        # simulations
        axes[row,col].plot(time_plot,sim_plot[:,-1],'b-',linewidth=1.5)
        if "WSE" in varname:
            station_lat = station_lat + [int(loc_lat)]
            station_lon = station_lon + [int(loc_lon)]
            station_idx  = station_idx  + [plot_idx]
        for time_ind in range(tstart_da,np.shape(wlv_da)[1]+tstart_da,1):
            x_time = np.arange(time_ind*dahour,time_ind*dahour+24,1)
            y_value = np.nanmean(wlv_da[:,int(time_ind-tstart_da),i,:],axis=0)
            plot_gradient_line(axes[row, col], x_time, y_value, 'Reds_r', linewidth=0.8,add_colorbar=(time_ind == tstart),fig=fig)        
        if rmse_da < rmse_sim:
            improve_rate[ind] = 1
        # axes[row,col].set_title(alabel[plot_idx],loc='left',fontsize=12)
        if col == 0:
            axes[row,col].set_ylabel(varname)
        if col == cols-1:
            ax2.set_ylabel('Rain (mm/h)')
        if row != 13:
            axes[row,col].set_xticklabels([],rotation=15)
        if ("WSE" in varname) & (row == 5):
            dates = pd.date_range(start='2024-07-24 15:00',end='2024-07-25 14:00',freq='6h')
            # axes[row,col].set_xticklabels(dates,rotation=25)
        plot_idx = plot_idx + 1
    for ind in range(plot_idx, rows_old * cols):
        row, col = divmod(ind, cols)
        fig.delaxes(axes[row, col])  # delete the empty 
    improve = np.nansum(improve_rate)/np.shape(np.where(~np.isnan(improve_rate)))[1]*100
    print('How many stations improve: ','{:.2f}'.format(improve),'%')
    if 'case' in pm_name:
        plt.savefig(exp_plotdir+'2.fore_case'+savename+'.jpg', format='jpg',dpi=600)
    else:
        plt.savefig(exp_plotdir+'2.fore_validation'+savename+'.jpg', format='jpg',dpi=600)
    plt.show()
    plt.close()
    return np.array(station_lat),np.array(station_lon),np.array(station_idx),None



# KGE, SVD, CC
def plot_foreline(ax,data_sim,data_da,ylabel):
    x = np.arange(np.shape(data_da)[0])
    if 'KGE' in ylabel:
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1.5)
    bar_width = 0.4
    # ax.axvline(x=3, color='gray', linestyle='--', linewidth=1.5)
    ax.plot(x,data_sim, color='blue', lw=2.2)
    ax.plot(x,data_da, color='red', lw=2.2)
    ax.set_xticks(np.arange(0,np.shape(data_da)[0]))
    if 'ils01' in pm_path or 'round' or 'all' in pm_path:
        ax.set_xticklabels([str(i) for i in range(1, 25)])
    if 'ils03' in pm_path:
        ax.set_xticklabels(['1-3','4-6','7-9','10-12','13-15','16-18','19-21','22-24'])
    if 'ils06' in pm_path:
        ax.set_xticklabels(['1-6','7-12','13-18','19-24'])
    if 'ils12' in pm_path:
        ax.set_xticklabels(['1-12','13-24'])
    ax.set_xlim(-0.5,np.nanmax(x)+0.5)
    ax.set_ylabel(ylabel,fontsize=12)
    ax.set_xlabel('Forecast Leading Time Before Peak Flow (h)',fontsize=12)
    
def plot_24lines(sim_matrix_wlv,da_matrix_wlv,varname):
    fig, axs = plt.subplots(3, 1, figsize=(8, 8),dpi=600)
    plot_foreline(axs[0],sim_matrix_wlv[:,5],da_matrix_wlv[:,5],'STD of ensembles') # standard daviation
    plot_foreline(axs[1],sim_matrix_wlv[:,4],da_matrix_wlv[:,:,4],'Pearson Correlation Coefficient')
    plot_foreline(axs[2],sim_matrix_wlv[:,2],da_matrix_wlv[:,:,2],'KGE')
    plt.tight_layout()
    plt.savefig(exp_plotdir+'5.24hour_line_'+varname+'.jpg', format='jpg',dpi=600)
    plt.show()

def plot_bar_error(ax,data_sim,data_da,ymin,ymax,ylabel,title):
    def get_stats(data):
        median = np.nanpercentile(data, 50, axis=0)
        q25 = np.nanpercentile(data, 25, axis=0)
        q75 = np.nanpercentile(data, 75, axis=0)
        q10 = np.nanpercentile(data, 10, axis=0)
        q90 = np.nanpercentile(data, 90, axis=0)
        err_lower = median - q10
        err_upper = q90 - median
        return median,q25,q75,err_lower,err_upper
    # sim
    med1 = np.nanpercentile(data_sim, 50)
    q25_1 = np.nanpercentile(data_sim, 25)
    q75_1 = np.nanpercentile(data_sim, 75)
    q10_1 = np.nanpercentile(data_sim, 10)
    q90_1 = np.nanpercentile(data_sim, 90)
    err1_low = med1 - q10_1
    err1_high = q90_1 - med1
    # DA
    med2,q25_2,q75_2, err2_low, err2_high = get_stats(data_da)
    x = np.arange(np.shape(data_da)[1]+1)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1.5)
    bar_width = 0.4
    ax.bar(x[0],q75_1-q25_1, bottom=q25_1, capsize=5, color='blue', width=bar_width, edgecolor='gray')
    ax.bar(x[1:], q75_2-q25_2, bottom=q25_2, capsize=5, color='red', width=bar_width, edgecolor='gray')
    for i in range(len(x)-1):
        ax.plot(x[i+1], med2[i], marker='_', color='black', markersize=12, linewidth=150, zorder=3)
    ax.plot(x[0], med1, marker='_', color='black', markersize=12, linewidth=150, zorder=3)
    ax.set_xticks(range(np.shape(data_da)[1]+1))
    if 'ils01' in pm_path:
        ax.set_xticklabels(['sim']+[str(i) for i in range(1, np.shape(data_da)[1]+1)])
    if 'ils03' in pm_path:
        ax.set_xticklabels(['sim','1-3','4-6','7-9','10-12','13-15','16-18','19-21','22-24'])
    if 'ils06' in pm_path:
        ax.set_xticklabels(['sim','1-6','7-12','13-18','19-24'])
    if 'ils12' in pm_path:
        ax.set_xticklabels(['sim','1-12','13-24'])
    ax.set_ylim(ymin,ymax)
    ax.set_xlim(-0.5,np.shape(data_da)[1]+1-0.5)
    ax.set_ylabel(ylabel,fontsize=12)
    ax.set_xlabel('Hours after DA (h)',fontsize=12)
    # ax.set_title('Case'+title,loc='left',fontsize=14)

def plot_24lines4(sim_matrix_wlv,da_matrix_wlv,varname,time_error_sim,time_error_da,exp_plotdir):
    fig, axs = plt.subplots(3, 1, figsize=(8, 8),dpi=600)
    plot_bar_error(axs[0],time_error_sim,time_error_da,-25,1,'Timing error (h)',' ') # timing errors
    plot_foreline(axs[1],sim_matrix_wlv[:,5],da_matrix_wlv[:,5],'STD of ensembles') # standard daviation
    plot_foreline(axs[2],sim_matrix_wlv[:,2],da_matrix_wlv[:,2],'KGE')
    legend_elements = [Patch(facecolor='blue', edgecolor='black', label='Open-loop experiment'),Patch(facecolor='red', edgecolor='black', label='DA experiment')]
    fig.legend(
        handles=legend_elements,
        loc='upper center',
        ncol=2,
        bbox_to_anchor=(0.5, 0.98),
        frameon=False)
    plt.tight_layout()
    plt.savefig(exp_plotdir+'5.24hour_4case2_'+varname+'.jpg', format='jpg',dpi=600)
    plt.show()
    if 'case2' in pm_name:
        np.save('sim_fore2024.npy',sim_matrix_wlv)
        np.save('da_fore2024.npy',da_matrix_wlv)
    if 'case1' in pm_name:
        np.save('sim_fore2022.npy',sim_matrix_wlv)
        np.save('da_fore2022.npy',da_matrix_wlv)

#%% draw the map
def draw_rmse(title,var_name,lon_plot,lat_plot,situ_rmse,vmin,vmax,lon,lat,exp_plotdir):
    fig = plt.figure(dpi = 600,figsize=(7.5,4.5))
    ax1 = fig.add_axes([0.1,0.1,0.8,0.8],projection=ccrs.PlateCarree())
    ax1.set_axis_off()
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
    gl.bottom_labels   = False
    gl.left_labels = False
    norm = mcolors.TwoSlopeNorm(vmin=min(situ_rmse), vcenter=0, vmax=-min(situ_rmse))
    if 'KGE' in title:
        bounds = [-1, -0.5, 0, 0.5, 0.75, 1]  
        ncolors = len(bounds) - 1
        cmap = plt.get_cmap('rainbow', ncolors)
        new_cmap = ListedColormap([cmap(i) for i in range(ncolors)])
        norm = BoundaryNorm(bounds, ncolors)
        con_plot2 = ax1.scatter(lon[lon_plot],lat[lat_plot],c=situ_rmse,s=7.5,zorder=2,cmap=new_cmap,norm=norm,transform=ccrs.PlateCarree(),alpha=0.8)
    else:
        con_plot2 = ax1.scatter(lon[lon_plot],lat[lat_plot],c=situ_rmse,s=7.5,zorder=2,cmap='coolwarm',transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,alpha=0.8)
    # station highlight  
    # ax1.scatter(lon[918],lat[556],c='black',alpha=0.5,zorder=3,transform=ccrs.PlateCarree(),marker='*',s=200)
    def color_bar(l,b,w,h):
      rect = [l,b,w,h]
      cbar_ax = fig.add_axes(rect)
      return cbar_ax
    [l1,b1,w1,h1] = [0.63,0.12,0.01,0.25]
    [l2,b2,w2,h2] = [0.63,0.32,0.01,0.25]
    cbar_ax2 = color_bar(l2,b2,w2,h2)
    if 'rivdph' in var_name:
        cb_title = '$(m)$'
    else:
        cb_title = ' '
    cb2 = plt.colorbar(con_plot2, cax=cbar_ax2,orientation="vertical",shrink=0.5)
    cb2.ax.set_title(cb_title,fontsize=8)
    ax1.set_title(title,loc='center')
    if 'round' in pm_path:
        plt.savefig(exp_plotdir+'3.map_'+title+pm_name[-1]+'.jpg', format='jpg',dpi=600)
    elif 'case' in pm_path:
        plt.savefig(exp_plotdir+'3.map_'+title+pm_name[-2:]+'.jpg', format='jpg',dpi=600)
    else:
        plt.savefig(exp_plotdir+'3.map_'+title+'.jpg', format='jpg',dpi=600)
    plt.show()
    plt.close()



def draw_each_figure(obs_lat,obs_lon,assim_select,id_info,varname,pm_path,inputdir,expname,dahour_str,pm_name,time_len,lon_select,lat_select,flood_time_select,data_obs_wlv,lon,lat,elvmean,simmean,fld_lonmin,fld_lonmax,fld_latmin,fld_latmax):
    situ_rmse  = np.full(np.shape(obs_lat)[0],np.nan)
    lat_plot   = np.full(np.shape(obs_lat)[0],1)
    lon_plot   = np.full(np.shape(obs_lat)[0],1)
    if varname == 'rivdph':
        typename = 'rivdph'
    else:
        typename = 'outflw'
    if 'all' in pm_path:
        if 'case' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+pm_name[7:]+'/'+typename+'/all/'
        else:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+typename+'/'
    else:
        if 'round' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/round'+pm_path[-4]+'/'+typename+'/'
        elif 'case' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+pm_name[7:]+'/'+typename+'/valid/'
        elif 'ens' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/MEPS/'+typename+'/'
        elif 'patch' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/patch/'+pm_name[7:]+'/'+typename+'/'
        else:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+typename+'/'
    # DA
    filemean = filepath+'meanA.bin'   # compare to the average of ensemble
    da_mean = read_bin_mean(filemean,time_len,lon_select,lat_select)  # compare with the open loop result
    # sim
    filemean = filepath+'meanC.bin'  # compare with average of ensembles
    # filemean = filepath+'oriC.bin' # compare with original simulations
    sim_mean = read_bin_mean(filemean,time_len,lon_select,lat_select)  # compare with the open loop result
    # if '00' in filename:
        # sim_mean = read_nc(filename,4)
    if varname == 'HQ':
        if 'round' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/round'+pm_path[-4]+'/rivdph/'
        else:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/rivdph/'
        da_dis_mean = read_bin_mean(filepath+'meanA.bin',time_len,lon_select,lat_select)
        Q_cal = cal_all_HQ(da_dis_mean,obs_select_wlv,obs_select_dis)
    for loc in range(0,np.shape(obs_lat)[0]):
        loc_lat = obs_lat[loc]
        loc_lon = obs_lon[loc]
        loc_assim = assim_select[loc]
        loc_id  = id_info[loc]
        station_flood_time = flood_time_select[loc,:]
        #obs data_obs_wlv/data_obs_dis
        if varname == 'rivdph':
            obs_grid = data_obs_wlv[:,loc_lat,loc_lon]
        else:
            obs_grid = data_obs_dis[:,loc_lat,loc_lon]
        if station_flood_time[0]<-900:
            # no flood event
            continue
        if (lon[loc_lon]<fld_lonmin) | (lon[loc_lon]>fld_lonmax) | (lat[loc_lat]<fld_latmin) | (lat[loc_lat]>fld_latmax):
            # exclude station out of flooding regions
            continue
        if np.all(obs_grid)<-900:
            continue
        else:
        # read data
            time1 = int(station_flood_time[0])
            time2 = int(station_flood_time[1]+1)
            if varname == 'rivdph':
                loc_obs = obs_grid[time1:time2]-elvmean[loc_lat,loc_lon]
                loc_sim = sim_mean[time1:time2,loc]-simmean[loc_lat,loc_lon]
                loc_da  = da_mean[time1:time2,loc]-simmean[loc_lat,loc_lon]
                obs_range = 0.01 # obs error range
            else:
                loc_obs = obs_grid[time1:time2]
                loc_sim = sim_mean[time1:time2,loc]
                obs_range = 10
                if varname == 'HQ':
                    loc_da  = Q_cal[time1:time2,loc]                
                else:
                    loc_da  = da_mean[time1:time2,loc]
            if np.all(loc_da[10:30]<10**(-8))==True:
                loc_da = np.full(len(da_mean[time1:time2,loc]),np.nan)
            # remove no obs stations 
            mse_sim = np.nanmean((loc_obs-loc_sim)**2)
            mse_da  = np.nanmean((loc_obs-loc_da)**2)
            if np.isnan(mse_sim) | np.isnan(mse_da) | np.all(np.isnan(loc_da)) | (np.nanmean(loc_obs)<obs_range):
                continue
            rmse = np.sqrt(mse_sim) - np.sqrt(mse_da)
            if (np.shape(loc_obs)[0]<25) | (np.shape(loc_sim)[0]<25) | (np.shape(loc_da)[0]<25):
                rmse = np.nan
                nrmse = np.nan
                imp = np.nan
            else:
                var_max  = np.nanmax(loc_obs)
                var_min  = np.nanmin(loc_obs)
                var_range = var_max-var_min
                if var_range ==0:
                    nrmse = np.nan
                else:
                    nrmse = rmse/var_range
                # exclude lake region
                if rmse<-0.:
                    if np.abs(elvmean[loc_lat,loc_lon]-simmean[loc_lat,loc_lon])>3:
                        rmse = np.nan
            # store nrmse, lat, loc
            situ_rmse[loc]=rmse
            lat_plot[loc]=int(loc_lat)
            lon_plot[loc]=int(loc_lon)
    return situ_rmse,lat_plot,lon_plot

def plot_eventbar(ax,sim1,sim2,sim3,da1,da2,da3,ylabel):
    x = np.arange(np.shape(sim1)[0])
    width=0.25
    ax.bar(x-width,sim1-da1,color='red',edgecolor='gray',width=width,label='2019 Flood Event')
    ax.bar(x,sim2-da2,color='blue',edgecolor='gray',width=width,label='2022 Flood Event')
    ax.bar(x+width,sim3-da3,color='orange',edgecolor='gray',width=width,label='2024 Flood Event')
    ax.set_xticks(np.arange(0,np.shape(sim1)[0]))
    ax.set_xlim(-0.5,np.nanmax(x)-0.5)
    if 'STE' in ylabel:
        ax.set_ylim(-1,5)
    ax.set_ylabel(ylabel,fontsize=12)
    ax.axhline(y=0,color='gray',linestyle='--')
    ax.set_xticks(np.arange(0,np.nanmax(x),1))
    ax.set_xticklabels(np.arange(1,np.nanmax(x)+1,1))
    ax.set_xlabel('Forecast Leading Time Before Peak Flow (h)',fontsize=12)

def plot_eventline(ax,sim1,sim2,sim3,da1,da2,da3,ylim,ylabel):
    x = np.arange(np.shape(sim1)[0])
    width=0.25
    ax.plot(x,sim1,color='red',linestyle='--',lw=1.25,label='2019 Flood OL')
    ax.plot(x,da1,color='red',lw=4.5,label='2019 Flood DA')
    ax.plot(x,sim2,color='blue',linestyle='--',lw=1.25,label='2022 Flood OL')
    ax.plot(x,da2,color='blue',lw=4.5,label='2022 Flood DA')
    ax.plot(x,sim3,color='orange',linestyle='--',lw=1.25,label='2024 Flood OL')
    ax.plot(x,da3,color='orange',lw=4.5,label='2024 Flood DA')
    ax.set_xticks(np.arange(0,np.shape(sim1)[0]))
    ax.set_xlim(-0.5,np.nanmax(x)-0.5)
    if 'STE' in ylabel:
        ax.set_ylim(-1,5)
    ax.set_ylabel(ylabel,fontsize=12)
    ax.axhline(y=0,color='gray',linestyle='--')
    ax.set_xticks(np.arange(0,np.nanmax(x),1))
    ax.set_xticklabels(np.arange(1,np.nanmax(x)+1,1))
    if 'STD' in ylabel:
        ax.set_ylim(0,ylim)
    else:
        ax.set_ylim(-ylim,0.2)
    ax.set_xlabel('Forecast Leading Time Before Peak Flow (h)',fontsize=12)
