import glob
import numpy as np
import os
import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_ils01 as pm    # out_msm_round*  # params_all_case3

def cal_assimilate_grid():
    with open(pm.CaMa_dir()+'/map/tej_01min/nextxy.bin', 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.int32)
        data_re  = np.reshape(data_bin,(2,1320,1500))
    flow_ix = data_re[0,:,:]-1  # next ix (ilon)
    flow_iy = data_re[1,:,:]-1  # next iy (ilat)
    flow_ix = np.where((flow_ix<0),np.nan,flow_ix)
    flow_iy = np.where((flow_iy<0),np.nan,flow_iy)
    river_grid_num = np.shape(np.where(~np.isnan(flow_ix)))[1]
    patchdir = '/data50/yingying/HydroDA/obs_patch/obs_patch06/obs_patch/'
    txt_files= glob.glob(patchdir+'*.txt')
    patch_num = np.shape(txt_files)[0]
    rate = patch_num/river_grid_num*100
    rateg = 1500/river_grid_num*100
    print('assimilated rate:','{:.2f}'.format(rate),'gauge rate:','{:.2f}'.format(rateg))
    patchname = [os.path.splitext(f)[0][-8:] for f in txt_files]
    with open(pm.CaMa_dir()+'/map/tej_01min/uparea.bin', 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(1320,1500))
    masked_data = np.full_like(data_re, np.nan)  
    valid = (~np.isnan(flow_ix)) & (~np.isnan(flow_iy))
    ix = flow_ix[valid].astype(int)
    iy = flow_iy[valid].astype(int)
    masked_data[iy, ix] = data_re[iy, ix]
    masked_data = np.where(masked_data<10**8,np.nan,masked_data)  # 100m2
    # [lat, lon, weight]
    data_weight = np.full((patch_num,3),np.nan)
    row = 0
    for ind in range(len(txt_files)):
        data  = np.array(np.loadtxt(txt_files[ind]))
        if data.ndim == 2:
            weight = np.nanmax(data[:,2])
        else:
            weight = data[2]
        lon_ind = int(patchname[ind][:4])-1
        lat_ind = int(patchname[ind][4:])-1
        if ~np.isnan(masked_data[lat_ind,lon_ind]):
            data_weight[row,1] = lon_ind 
            data_weight[row,0] = lat_ind 
            data_weight[row,2] = weight
            row = row+1
        
    return data_weight
data_weight = cal_assimilate_grid()
data_new = data_weight[~np.all(np.isnan(data_weight),axis=1)]
np.save('weight.npy',data_new)
