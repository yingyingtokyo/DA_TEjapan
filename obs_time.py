import numpy as np
import os
import re
from datetime import datetime, timedelta

import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_all_case1 as pm    # out_msm_round*  # params_all_case3
import importlib
importlib.reload(pm)

stime = datetime(2022,1,1,0)
etime = datetime(2022,12,31,23)
time_len = int((etime - stime).total_seconds() / 3600) + 1
def read_bin2d(fname):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(1320,1500))
    return data_re

def write_bin2d(fname, result):
    roff32 = result.astype(np.float32)
    roff_flatten = roff32.flatten()
    with open(fname, 'wb') as f:
        roff_flatten.tofile(f)

def read_obs(yyyy):
    current_time = stime
    data_obs = np.full((time_len,1320,1500),np.nan)
    ind = 0
    while current_time<etime:
        ctime = current_time.strftime('%Y%m%d%H')
        fname = '/data42/yingying/obs/'+yyyy+'UST/wlv/'+ctime+'.bin'
        data_obs[ind,:,:] = read_bin2d(fname)
        current_time += timedelta(hours=1)
        ind = ind + 1
    data_obs = np.where(data_obs<-900,np.nan,data_obs)
    return data_obs

data_obs = read_obs('2022')
# [time_len,1320,1500]
data_obs_new = np.full((np.shape(data_obs)[0],1320,1500),np.nan)
# 9 hours time shift (but only for 2022)
data_obs_new[:-9,:,:] = data_obs[9:,:,:]

def write_obs(yyyy):
    current_time = stime
    ind = 0
    while current_time<etime:
        ctime = current_time.strftime('%Y%m%d%H')
        fname = '/data42/yingying/obs/'+yyyy+'UST_new/wlv/'+ctime+'.bin'
        write_bin2d(fname,data_obs_new[ind,:,:])
        current_time += timedelta(hours=1)
        ind = ind + 1
write_obs('2022')
