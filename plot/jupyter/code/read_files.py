import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case21 as pm
pm_path = pm.__file__
pm_name = pm.__name__

import numpy as np
import calendar
import os

import importlib
import basic_command as basic
importlib.reload(basic)
from basic_command import find_id, read_bin2d, read_bin3d, read_bin_int, read_bin_station, read_nc

def all_select_plot(pm_name):
    savepath = '/data50/yingying/HydroDA/plot/exp_ils01/' + pm_name[-6:-1]
    lat_all_flood = []
    lon_all_flood = []
    assim_all_flood = []
    obs_station_all = []
    da_station_all  = []
    sim_station_all = []
    for case_ind in range(1, 6):
        lat_each_flood = np.load(savepath + f'{case_ind}lat.npy')        # shape: (24, 20)
        lon_each_flood = np.load(savepath + f'{case_ind}lon.npy')        # shape: (24, 20)
        assim_num_flood = np.load(savepath + f'{case_ind}assim_num.npy') # shape: (24, 20)
        obs_flood = np.load(savepath + f'{case_ind}obs_station.npy')     # shape: (24, 8)
        da_flood  = np.load(savepath + f'{case_ind}da_station.npy')      # shape: (24, 20, 8)
        sim_flood = np.load(savepath + f'{case_ind}sim_station.npy')     # shape: (24, 20, 8)
        lat_all_flood.append(lat_each_flood)
        lon_all_flood.append(lon_each_flood)
        assim_all_flood.append(assim_num_flood)
        obs_station_all.append(obs_flood)
        da_station_all.append(da_flood)
        sim_station_all.append(sim_flood)
    lat_all_flood      = np.concatenate(lat_all_flood, axis=0)  # shape: (5×24, 20)
    lon_all_flood      = np.concatenate(lon_all_flood, axis=0)
    assim_all_flood    = np.concatenate(assim_all_flood, axis=0)
    obs_station_all    = np.concatenate(obs_station_all, axis=-1)  # shape: (24, 20, 20)
    da_station_all     = np.concatenate(da_station_all, axis=-1)
    sim_station_all    = np.concatenate(sim_station_all, axis=-1)
    return lat_all_flood, lon_all_flood, assim_all_flood, obs_station_all, da_station_all, sim_station_all

def def_oscillating(obs_wlv,time_len):
    min_crossings=time_len/8
    obs_mean = np.nanmean(obs_wlv)
    anomaly  = obs_wlv-obs_mean
    signs    = np.sign(anomaly)
    crossings = np.sum(signs[1:] * signs[:-1]<0)
    # if True: remove this observation
    return crossings >= min_crossings
    
def def_flood_time(year,ix,iy,id_info,lon_point_ind,lat_point_ind,time_len,data_obs_wlv,elvmean,tstart):
    # [start_time, end_time]
    flood_time = np.full((len(lon_point_ind),2),-999)
    for sta_ind in range(0,len(lon_point_ind)):
        sta_id = find_id(lat_point_ind[sta_ind],lon_point_ind[sta_ind],ix,iy,id_info)
        filen = str(sta_id)+'.bin'
        delta_all = np.array([])
        # station obs
        obs_flood_wlv = data_obs_wlv[:,lat_point_ind[sta_ind],lon_point_ind[sta_ind]]
        delta_obs = np.abs(obs_flood_wlv[1:]-obs_flood_wlv[:-1])
        delta_obs = np.append(delta_obs,np.nan)
        if_osc= def_oscillating(obs_flood_wlv,time_len)
        if if_osc == False:
            # observation is not oscillating
            # annual obs at stations
            filepath = pm.DA_dir()+'/Empirical_LocalPatch/MLIT/data/data_wlv_hour'+"%04d"%(year)+'/'
            if os.path.exists(filepath+filen):
                with open(filepath+filen, 'rb') as file:
                    data_bin = np.fromfile(file, dtype=np.float64)
                data_bin = np.where((data_bin>10**4)|(data_bin<-900),np.nan,data_bin)
                if calendar.isleap(year):
                    days = 366
                else:
                    days = 365
                data_re  = np.reshape(data_bin,days*24)
                delta_data = data_re[1:]-data_re[:-1]
                delta_all   = np.concatenate((delta_all,delta_data))
                delta = delta_all[~np.isnan(delta_all)]
                if np.shape(delta)[0]>100:
                    thold_value = np.percentile(delta,99)
                    mean_value  = np.nanmean(delta)-elvmean[lat_point_ind[sta_ind],lon_point_ind[sta_ind]]
                else:
                    thold_value = 10**5
                    mean_value  = 10**5                
            else:
                thold_value = 10**5
                mean_value  = 10**5
            ftime = np.where((delta_obs>thold_value)&(obs_flood_wlv>mean_value))[0]
            ftime = ftime[ftime>tstart]
            if (np.shape(ftime)[0]>1):
                # flood_occur
                obsmax_ind = np.nanargmax(obs_flood_wlv)
                if (obsmax_ind<ftime[0]):
                    ftime_mean = np.where(obs_flood_wlv>mean_value)[0]
                    flood_time[sta_ind,0] = int(ftime_mean[0])
                    flood_time[sta_ind,1] = int(ftime[-1])            
                elif (obsmax_ind>ftime[-1]):
                    ftime_mean = np.where(obs_flood_wlv>mean_value)[0]
                    flood_time[sta_ind,0] = int(ftime[0])
                    flood_time[sta_ind,1] = int(ftime_mean[-1])                       
                else:
                    flood_time[sta_ind,0] = int(ftime[0])
                    flood_time[sta_ind,1] = int(ftime[-1])            
    return flood_time

def cal_all_HQ(da_station_wlv,obs_station,obs_station_dis):
    cal = 0
    Q_obs = np.full((np.shape(da_station_wlv)[0],np.shape(da_station_wlv)[1]),np.nan)
    for station_ind in range(0,np.shape(obs_station)[1]):
        if (np.shape(np.unique(obs_station_dis[tstart:,station_ind]))[0]<((time_len-tstart)/4)) | (np.shape(np.unique(obs_station[tstart:,station_ind]))[0]<((time_len-tstart)/4)):
            continue
        obs_wlv = obs_station[tstart:,station_ind]
        obs_dis = obs_station_dis[tstart:,station_ind]
        mask = ~np.isnan(obs_wlv) & ~np.isnan(obs_dis)
        obs_wlv_new = obs_wlv[mask]
        obs_dis_new = obs_dis[mask]
        coeffi = np.polyfit(obs_wlv_new,obs_dis_new,2) # second-order polynomic
        a, b, c = coeffi
        Q_fit  = a*obs_wlv_new**2 + b*obs_wlv_new + c
        R,_ = pearsonr(obs_dis_new,Q_fit)
        da_wlv  = da_station_wlv[tstart:,station_ind]
        Q_obs[tstart:,station_ind]  = a*da_wlv**2 + b*da_wlv + c
        if R<0.85:
            continue
        cal = cal + 1
    return Q_obs

def obs_select(obs_lat,obs_lon,data_obs_wlv):
    obs_wlv = np.full((np.shape(data_obs_wlv)[0],np.shape(obs_lat)[0]),np.nan)
    for loc in range(0,np.shape(obs_lat)[0]):
        loc_lat = obs_lat[loc]
        loc_lon = obs_lon[loc]
        obs_wlv[:,loc] = data_obs_wlv[:,loc_lat,loc_lon]
    return obs_wlv

def performance(situ_rmse):
    station_num = np.shape(np.where(situ_rmse>0))[1]+np.shape(np.where(situ_rmse<0))[1]
    if (station_num) == 0:
        nmse_pos = np.nan
    else:
        nmse_pos = np.shape(np.where(situ_rmse>0))[1]/station_num
    nmse_neg = 1-nmse_pos
    print('station numbers:',station_num)
    print('better',nmse_pos*100,'worse',nmse_neg*100)
    print('-------')


def performance_all(situ_rmse,situ_nrmse,situ_imp,situ_nse,situ_kge):
    station_num = np.shape(np.where(situ_nrmse>0))[1]+np.shape(np.where(situ_nrmse<0))[1]
    if (station_num) == 0:
        nrmse_pos = np.nan
    else:
        nrmse_pos = np.shape(np.where(situ_nrmse>0))[1]/station_num
    nrmse_neg = 1-nrmse_pos
    nse_da_good  = np.shape(np.where(situ_nse[:,0]>0.5))[1]
    nse_sim_good = np.shape(np.where(situ_nse[:,1]>0.5))[1]
    kge_da_good  = np.shape(np.where(situ_kge[:,0]>0.25))[1]
    kge_sim_good = np.shape(np.where(situ_kge[:,1]>0.25))[1]
    print('station numbers:',station_num)
    print('better',nrmse_pos*100,'worse',nrmse_neg*100)
    print('average performance',np.nanmean(situ_rmse),np.nanmean(situ_nrmse),np.nanmean(situ_imp)*100)
    print('nse_da_good: ',nse_da_good/station_num*100,'nse_sim_good',nse_sim_good/station_num*100)
    print('kge_da_good: ',kge_da_good/station_num*100,'kge_sim_good',kge_sim_good/station_num*100)
    print('-------')
    values = [nrmse_pos*100, nrmse_neg*100, np.nanmean(situ_rmse), np.nanmean(situ_nrmse), np.nanmean(situ_imp)*100]

