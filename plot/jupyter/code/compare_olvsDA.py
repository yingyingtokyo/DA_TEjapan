import sys
params_path = '/data50/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_case21 as pm
pm_path = pm.__file__
pm_name = pm.__name__

import numpy as np
import netCDF4 as nc

import importlib
import basic_command as basic
importlib.reload(basic)
from basic_command import find_id, read_bin2d, read_bin3d, read_bin_int,read_bin_station


def read_loc(file_valid_lat,file_valid_lon,obs_situ_dir,num_assim_path,latx_valid,lony_valid,assim_num,lat_select_file,lon_select_file):
    # Load validation index files
    latx_valid = np.loadtxt(obs_situ_dir + file_valid_lat, dtype=int)
    lony_valid = np.loadtxt(obs_situ_dir + file_valid_lon, dtype=int)
    # Load assimilation count
    assim_pixel = read_bin_int(num_assim_path, 2)
    assim_plot_num = assim_pixel[0, :, :]
    assim_num = assim_plot_num[latx_valid - 1, lony_valid - 1]
    # Filter out zero-assimilation points
    mask = assim_num != 0
    lat_select_file = np.delete(latx_valid,assim_num==0)
    lon_select_file = np.delete(lony_valid,assim_num==0)
    # Remove duplicates
    coords = np.stack((lat_select_file,lon_select_file),axis=1)
    _,unique_indices = np.unique(coords,axis=0,return_index = True)  # delete repeated elements
    unique_ind_sorted = np.sort(unique_indices)
    lat_unique = lat_select_file[unique_ind_sorted]
    lon_unique = lon_select_file[unique_ind_sorted]
    sorted_ind = np.lexsort((lon_unique, lat_unique))
    lat_valid_plot = lat_unique[sorted_ind]-1
    lon_valid_plot = lon_unique[sorted_ind]-1
    return sorted_ind,lat_valid_plot,lon_valid_plot


# read results from all rounds
def read_all_round(varnames,savename,path_var,time_len,expname,dahour_str,inputdir,lat_point_ind,lon_point_ind):
    rmse_all = np.empty((0, 8))
    da_station  = np.full((time_len,pm.ens_mem(),500),np.nan)
    sim_station = np.full((time_len,pm.ens_mem(),500),np.nan)
    off_station = np.full((time_len,500),np.nan)
    count = 0
    count_row = 0
    for var in varnames:
        # [plot_idx, rmse_da, rmse_sim, kge_da, kge_sim, lat, lon, assim_num]
        if 'round' in var:
            rmse_each = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse_round'+savename+var[-1]+'.npy')
            filepath2 = 'valid/round'+var[-1]+'/'
            path_off = inputdir+'/exp_'+expname+dahour_str+'/round'+var[-1]+'/rivdph/'
            lat_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lat_point_round'+var[-1]+'.npy')
            lon_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lon_point_round'+var[-1]+'.npy')
        elif 'case' in var:
            rmse_each = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse_case'+savename+var[-1]+'.npy')
            filepath2 = var[-1]+'/valid/'
            path_off = inputdir+'/exp_'+expname+dahour_str+'/round'+var[-1]+'/rivdph/'
            lat_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lat_case_ind'+var[-2]+var[-1]+'.npy')
            lon_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lon_case_ind'+var[-2]+var[-1]+'.npy')
        else:
            rmse_each = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse'+savename+'.npy')
            filepath2 = 'valid/'
            path_off = inputdir+'/exp_'+expname+dahour_str+'/rivdph/'
            lat_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lat_point_ind.npy')
            lon_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lon_point_ind.npy')
        
        # offline DA 
        fileoff = path_off+'offC.bin'
        with open(fileoff, 'rb') as file:
            data_bin = np.fromfile(file, dtype=np.float32)
            data_off  = np.reshape(data_bin,(time_len,len(lat_point_ind)))

        # row_ind = list()
        for sta_ind in range(len(lat_point_ind)):
        #     row_ind = row_ind + [int(sorted_ind[sta_ind])]
        # for sta_ind in ind_keep:
            # print(row_ind[sta_ind])
            off_station[:,count_row] = data_off[:,sta_ind]
            count_row = count_row + 1
            
        n = int(len(rmse_each) / 8)
        rmse_each = rmse_each.reshape(n, 8)
        rmse_all  = np.vstack([rmse_all,rmse_each])
        filepath = path_var+filepath2  # riv_outdir + filepath2 or dis_outdir + filepath2
        for sta_ind in range(np.shape(rmse_each)[0]):
            lat_point = int(rmse_each[sta_ind,5])
            lon_point = int(rmse_each[sta_ind,6])
            mean_station = simmean[lat_point,lon_point]
            for ens in range(0,pm.ens_mem()):
                fname = '{:04d}'.format(lat_point+1) +  '{:04d}'.format(lon_point+1) + '{:02d}'.format(ens+1)
                da_each = read_bin_station(filepath+fname + 'A.bin')
                sim_each= read_bin_station(filepath+fname + 'C.bin')
                if 'wlv' in savename:
                    da_station[:,ens,count] = da_each-mean_station
                    sim_station[:,ens,count]= sim_each-mean_station
                else:
                    da_station[:,ens,count] = da_each
                    sim_station[:,ens,count]= sim_each
            count = count + 1
    result_da  = da_station[:,:,:count]
    result_sim = sim_station[:,:,:count]
    result_off = off_station[:,:count_row]
    return result_da, result_sim, result_off, np.array(rmse_all),
