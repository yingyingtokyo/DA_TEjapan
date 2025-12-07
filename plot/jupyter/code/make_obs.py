import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case21 as pm
pm_path = pm.__file__
pm_name = pm.__name__

import numpy as np
from datetime import datetime, timedelta

import basic_command as basic
from basic_command import find_id, read_bin2d, read_bin3d, read_bin_int, read_bin_station, read_nc


def read_obs(var,stime,etime,time_len,obsdir):
    current_time = stime
    data_obs = np.full((time_len,1320,1500),np.nan)
    ind = 0
    while current_time<etime:
        ctime = current_time.strftime('%Y%m%d%H')
        if var=='wlv':
            fname = obsdir+'/'+ctime+'.bin'
        else:
            fname = '/data42/yingying/obs/2019UST/dis/'+ctime+'.bin'
        data_obs[ind,:,:] = read_bin2d(fname)
        current_time += timedelta(hours=1)
        ind = ind + 1
    data_obs = np.where(data_obs<-900,np.nan,data_obs)
    return data_obs

def read_stations(var,pm_path,stime,etime,time_len,obsdir,elvmean,inputdir,expname,dahour_str,lat_valid_plot,lon_valid_plot,rainf,simmean):
    if 'all' in pm_path:
        # assimilated
        if 'case' in pm_path:
            filepath2 = pm_name[7:12]
        else:
            filepath2 = 'all/'
    else:
        # validated
        if 'round' in pm_path:
            filepath2 = 'valid/round'+pm_path[-4]+'/'
        elif 'case' in pm_path:
            filepath2 = pm_name[7:12]
        elif 'ens' in pm_path:
            filepath2 = 'ens/'
        else:
            filepath2 = 'valid/'
    if var == 'rivdph':
        data_obs_wlv = read_obs('wlv',stime,etime,time_len,obsdir)
        data_obs = data_obs_wlv - elvmean
        if 'case' in pm_path:
            if 'all' in pm_path:
                filepath = inputdir + '/exp_' + expname + dahour_str + '/rivdph/' + filepath2 + 'all/'
            else:
                filepath = inputdir + '/exp_' + expname + dahour_str + '/rivdph/' + filepath2 + '/valid'+str(pm.val_round())+'/'
        else:
            filepath = inputdir + filepath2
    else:
        data_obs = read_obs('dis',stime,etime,time_len,obsdir)
        if 'case' in pm_path:
            if 'all' in pm_path:
                filepath = inputdir + '/exp_' + expname + dahour_str + '/outflw/' + filepath2 + 'all/'
            else:
                filepath = inputdir + '/exp_' + expname + dahour_str + '/outflw/' + filepath2 + '/valid'+str(pm.val_round())+'/'
        else:
            filepath = inputdir + filepath2
    da_station = np.full((time_len,pm.ens_mem(),len(lon_valid_plot)),np.nan)
    sim_station = np.full((time_len,pm.ens_mem(),len(lon_valid_plot)),np.nan)
    obs_station = np.full((time_len,len(lon_valid_plot)),np.nan)
    mean_station = np.full((len(lon_valid_plot)),np.nan)
    prep_station = np.full((time_len,len(lon_valid_plot)),np.nan)
    for sta_ind in range(len(lon_valid_plot)):
        obs_station[:,sta_ind]  = data_obs[:,int(lat_valid_plot[sta_ind])-1,int(lon_valid_plot[sta_ind])-1]
        prep_station[:,sta_ind] = rainf[:,int(lat_valid_plot[sta_ind])-1,int(lon_valid_plot[sta_ind])-1]
        mean_station[sta_ind] = simmean[int(lat_valid_plot[sta_ind])-1,int(lon_valid_plot[sta_ind])-1]
        for ens in range(0,pm.ens_mem()):
            fname = '{:04d}'.format(int(lat_valid_plot[sta_ind])) +  '{:04d}'.format(int(lon_valid_plot[sta_ind])) + '{:02d}'.format(ens+1)
            da_station[:,ens,sta_ind] = read_bin_station(filepath+fname + 'A.bin')
            sim_station[:,ens,sta_ind] = read_bin_station(filepath+fname + 'C.bin')
    return data_obs,da_station,sim_station,obs_station,mean_station,prep_station
