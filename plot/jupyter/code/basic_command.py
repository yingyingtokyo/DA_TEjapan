import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case21 as pm
pm_path = pm.__file__
pm_name = pm.__name__

import netCDF4 as nc
import numpy as np
import glob
import os
import rasterio
import matplotlib.ticker as mticker
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from scipy.signal import correlate

def find_id(situ_lat,situ_lon,ix,iy,id_info):
    tol = 1e-6
    matched_row = np.where((np.abs(ix-situ_lon)<tol) & (np.abs(iy-situ_lat)<tol))[0][0]
    station_id  = id_info[matched_row]
    return station_id
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
def read_prep(file):
    nf = nc.Dataset(file,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    print(varname)
    lon = np.array(nf.variables[varname[1]][:])
    lat = np.array(nf.variables[varname[2]][:])
    var = np.array(nf.variables[varname[3]][:])
    var = var*3600 # unit: mm/h 
    var = np.array(np.where(var<-900,np.nan,var))
    return lat,lon,var
def read_bin_mean(fname,time_len,lon_select,lat_select):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(time_len,len(lat_select)))
    return data_re
def read_bin_assim(fname,time_len_assim,lon_select,lat_select):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(time_len_assim,len(lat_select)))
    return data_re

def cal_rmse(y_true,y_pred,varname,time_len):
    len_obs = np.shape(np.where(~np.isnan(y_true)))[1]
    if "Outflow" in varname:
        if np.shape(np.where(~np.isnan(y_true[24*2:24*4])))[1]<10:
            rmse = np.nan
        else:
            rmse = np.sqrt(np.nanmean((y_pred-y_true)**2))
    if "WSE" in varname:
        if len_obs<int(time_len)/4:
            rmse = np.nan
        else:
            rmse = np.sqrt(np.nanmean((y_pred-y_true)**2))
    return rmse
def cal_kge(sim,obs,varname):
    sim_nonan = sim[(~np.isnan(sim))&(~np.isnan(obs))]
    obs_nonan = obs[(~np.isnan(sim))&(~np.isnan(obs))]
    r = np.corrcoef(sim_nonan, obs_nonan)[0, 1]
    alpha = np.std(sim_nonan) / np.std(obs_nonan)
    beta = np.mean(sim_nonan) / np.mean(obs_nonan)
    kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    return kge
def cal_nse(sim,obs,varname):
    sim_nonan = sim[(~np.isnan(sim))&(~np.isnan(obs))]
    obs_nonan = obs[(~np.isnan(sim))&(~np.isnan(obs))]
    numerator = np.sum((obs_nonan - sim_nonan) ** 2)
    denominator = np.sum((obs_nonan - np.mean(obs_nonan)) ** 2)
    return 1 - (numerator / denominator)
def cal_single_metric(obs_data,sim_data):
    if np.shape(obs_data)[0] != np.shape(sim_data)[0]:
        print('shape:',np.shape(obs_data)[0],np.shape(sim_data)[0])
    mask = ~np.isnan(obs_data) & ~np.isnan(sim_data)
    if np.shape(sim_data)[0]>20:
        obs = obs_data[mask]
        sim = sim_data[mask]
    else:
        obs = obs_data
        sim = sim_data
    if ((np.all(mask) == True) & (np.shape(sim)[0]>1)) | (np.shape(sim)[0]>20):
        rmse = np.sqrt(np.mean((sim-obs)**2))
        bias = np.mean(sim-obs)
        # kge
        r = np.corrcoef(sim, obs)[0, 1]
        alpha = np.std(sim) / np.std(obs)
        beta = np.mean(sim) / np.mean(obs)
        kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
        # r
        pearson_r, _ = pearsonr(obs, sim)
        # Timing error (lag of peak)
        numerator = np.sum((obs - sim) ** 2)
        denominator = np.sum((obs - np.mean(obs)) ** 2)
        # nse
        nse = 1 - (numerator / denominator)
        return rmse,bias,kge,nse,pearson_r
    else:
        return np.nan,np.nan,np.nan,np.nan,np.nan
def cal_peak_matric(obs,sim):
    rmse = np.sqrt(np.mean((sim-obs)**2))
    bias = np.mean(sim-obs)
    # kge
    r = np.corrcoef(sim, obs)[0, 1]
    alpha = np.std(sim) / np.std(obs)
    beta = np.mean(sim) / np.mean(obs)
    kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    # r
    pearson_r, _ = pearsonr(obs, sim)
    # Timing error (lag of peak)
    numerator = np.sum((obs - sim) ** 2)
    denominator = np.sum((obs - np.mean(obs)) ** 2)
    # nse
    nse = 1 - (numerator / denominator)
    return rmse,bias,kge,nse,pearson_r


def evaluate_timeerror(wlv_da,dahour,obs_point_ind,sim_point_ind,tstart_fore,obs_station_all,lat_all_flood,lon_all_flood,simmean,sim_station_all,da_station_all,time_len):
    time_date = np.arange(0,23,dahour)
    time_error_all = np.full((np.shape(obs_point_ind)[0],len(time_date)),np.nan) # from 1 to 24 hour leading time
    time_error_sim = np.full((np.shape(obs_point_ind)[0]),np.nan) # from 1 to 24 hour leading time
    plot_idx = 0
    for station_ind in range(0,np.shape(sim_point_ind)[0]):
        obs_data = obs_station_all[tstart_fore:,station_ind]
        loc_lat  = lat_all_flood[station_ind]-1
        loc_lon  = lon_all_flood[station_ind]-1
        mean_station = simmean[int(loc_lat),int(loc_lon)]
        sim_plot = sim_station_all[tstart_fore:,:,station_ind]-mean_station
        da_plot  = da_station_all[tstart_fore:,:,station_ind]-mean_station
        # RMSE
        rmse_da = cal_rmse(obs_data,np.nanmean(da_plot,axis=1),'WSE',time_len)
        rmse_sim= cal_rmse(obs_data,np.nanmean(sim_plot,axis=1),'WSE',time_len)
        if (np.isnan(rmse_da)) | (np.isnan(rmse_sim)):
            continue
        if (np.shape(np.unique(da_plot))[0]<2):  # exclude the fail DA result
            continue
        if (np.abs(np.nanmax(obs_data))<0.1):  # exclude the gauges at small rivers
            continue
        obs_peak = np.argmax(obs_data) + tstart_fore
        # original simulation
        sim_data = sim_point_ind[plot_idx,:,-1]
        sim_peak = np.argmax(sim_data) + tstart_fore
        time_error_sim[plot_idx]  = sim_peak-obs_peak
        obs_peak_ind = int((obs_peak-1)/dahour)
        obs_ind_remain = obs_peak-obs_peak_ind*dahour
        da_data = np.nanmean(wlv_da[:,:,station_ind,:],axis=0)
        for time_ind in range(0,23,dahour):
            sim_peak = np.argmax(da_data[int(time_ind/dahour),:])
            time_error_all[plot_idx,int(time_ind/dahour)] = sim_peak-(24-int(time_ind/dahour)*dahour-obs_ind_remain)#-tstart_fore
            # <0 earlier then real
            # >0 later then real
        plot_idx = plot_idx+1
    # [stations, leadtime]
    time_error_mean = np.nanmean(time_error_all,axis=0)
    return time_error_all[:,::-1],time_error_sim
#%% evaluate the 24 hours forecast performance which is ahead of peak flow
def cal_STD(sim_ens):
    # sim_ens = np.nan_to_num(sim_ens, nan=0.0, posinf=0.0, neginf=0.0)
    std_time = np.std(sim_ens)
    return std_time
def evaluate_hour(obs_station,sim_station,wlv_da,hour_str,time_leng,lat_all_flood,lon_all_flood,obs_station_all,da_station_all,sim_station_all,tstart_fore,simmean,time_len,wlv_sim):
    valid_indices = [i for i in range(len(lat_all_flood))]
    matrix_da  = np.full((len(lat_all_flood),5),np.nan)
    matrix_sim = np.full((len(lat_all_flood),5),np.nan)
    matrix_da_res  = np.full(7,np.nan)
    matrix_sim_res = np.full(7,np.nan)
    # plot stations
    plot_idx = 0
    for ind, i in enumerate(valid_indices):
        obs_plot = obs_station_all[tstart_fore:,i]
        loc_lat  = lat_all_flood[i]-1
        loc_lon  = lon_all_flood[i]-1
        mean_station = simmean[int(loc_lat),int(loc_lon)]
        sim_plot = sim_station_all[tstart_fore:,:,i]-mean_station
        da_plot  = da_station_all[tstart_fore:,:,i]-mean_station
        # RMSE
        rmse_da = cal_rmse(obs_plot,np.nanmean(da_plot,axis=1),'WSE',time_len)
        rmse_sim= cal_rmse(obs_plot,sim_plot[:,-1],'WSE',time_len)  
        if (np.isnan(rmse_da)) | (np.isnan(rmse_sim)):
            continue
        if (np.shape(np.unique(da_plot))[0]<2):  # exclude the fail DA result
            continue
        if (np.abs(np.nanmax(obs_plot))<0.1):  # exclude the gauges at small rivers
            continue  
        # da_plot  = np.nanmean(wlv_da[:,:,plot_idx,int(hour_str)-25],axis=0)
        # sim_plot = np.nanmean(wlv_sim[:,:,plot_idx,int(hour_str)-25],axis=0)
        da_plot  = wlv_da[0,:,plot_idx,int(hour_str)-25]
        sim_plot = wlv_sim[0,:,plot_idx,int(hour_str)-25]
        da_sta   = wlv_da[:,:,plot_idx,int(hour_str)-25]
        sim_sta  = wlv_sim[:,:,plot_idx,int(hour_str)-25]
        time_ind = np.argmax(obs_plot)   # peak flow time
        time_sind= time_ind - time_leng
        time_eind= time_ind + 2
        if time_ind<2:  # fake flood event
            continue
        obs_cal  = obs_plot[time_sind:time_eind]
        ######## only the 1h lead-time result has been considered  ##########
        # sim_cal  = sim_plot[time_sind:time_eind,-1] 
        # matrix_sim[5]  = cal_STD(sim_plot[time_sind:time_eind,:])
        #####################################################################
        # wlv [ens,time,sta_ind,lead_time]
        # important to consider the length of forecast time
        da_cal   = da_plot[time_sind-int(hour_str):time_eind-int(hour_str)]
        sim_cal  = sim_plot[time_sind-int(hour_str):time_eind-int(hour_str)]
        if (np.shape(obs_cal)[0]<3) | (np.shape(sim_cal)[0]<3) | (np.shape(da_cal)[0]<3):
            continue
        if (np.shape(obs_cal)[0]!=np.shape(sim_cal)[0]) | (np.shape(da_cal)[0]!=np.shape(obs_cal)[0]) | (np.shape(sim_cal)[0]!=np.shape(da_cal)[0]):
            continue
        matrix_da[plot_idx,0],matrix_da[plot_idx,1],matrix_da[plot_idx,2],matrix_da[plot_idx,3],matrix_da[plot_idx,4] = cal_peak_matric(obs_cal,da_cal)
        matrix_sim[plot_idx,0],matrix_sim[plot_idx,1],matrix_sim[plot_idx,2],matrix_sim[plot_idx,3],matrix_sim[plot_idx,4] = cal_peak_matric(obs_cal,sim_cal)
        # calculate STD
        matrix_da_res[5]   = cal_STD(da_sta[:,time_sind-int(hour_str):time_eind-int(hour_str)])
        matrix_sim_res[5]  = cal_STD(sim_sta[:,time_sind-int(hour_str):time_eind-int(hour_str)])
        plot_idx = plot_idx + 1
    matrix_da_res[:5] =np.nanmean(matrix_da,axis=0)
    matrix_sim_res[:5]=np.nanmean(matrix_sim,axis=0)
    # calculate how many stations with KGE>0.5
    matrix_da_res[6]  = np.shape(np.where(matrix_da[:,2]>0.5))[1]
    matrix_sim_res[6] = np.shape(np.where(matrix_sim[:,2]>0.5))[1]
    return matrix_da_res, matrix_sim_res
def evaluate_tot_hour(obs_station,sim_station,wlv_da,lat_all_flood,lon_all_flood,obs_station_all,da_station_all,sim_station_all,tstart_fore,simmean,time_len,wlv_sim):
    # matrix [rmse,bias,kge,nse,pearson_r,STD]
    da_matrix_wlv  = np.full((24,7),np.nan)
    sim_matrix_wlv = np.full((24,7),np.nan)
    time_leng = 1
    for hour_ind in range(1,1+24):
        hour_str = str(hour_ind)
        da_matrix_wlv[hour_ind-1,:], sim_matrix_wlv[hour_ind-1,:] = evaluate_hour(obs_station,sim_station,wlv_da,hour_str,time_leng,lat_all_flood,lon_all_flood,obs_station_all,da_station_all,sim_station_all,tstart_fore,simmean,time_len,wlv_sim)   
    # print(sim_matrix_wlv[0,2])
    return da_matrix_wlv, sim_matrix_wlv
    
def save_select_points(pm_name,lat_valid_plot,lon_valid_plot,assim_num_plot):
    np.save('/data50/yingying/HydroDA/plot/exp_ils01/'+pm_name[-6:]+'lat.npy',lat_valid_plot)
    np.save('/data50/yingying/HydroDA/plot/exp_ils01/'+pm_name[-6:]+'lon.npy',lon_valid_plot)
    np.save('/data50/yingying/HydroDA/plot/exp_ils01/'+pm_name[-6:]+'assim_num.npy',assim_num_plot)

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
    
def low_resolution(var,size,lat_min,lat_max,lon_min,lon_max):
    var_low = np.full((int((lat_max-lat_min)*3600/size),int((lon_max-lon_min)*3600/size)),np.nan)
    for row in range(0,np.shape(var_low)[0]):
        if (np.all(np.isnan(var[row*size:(row+1)*size,:]))==1):
            continue
        for col in range(0,np.shape(var_low)[1]):
            var_low[row,col] = np.nanmean(var[row*size:(row+1)*size,col*size:(col+1)*size])
    return var_low 
