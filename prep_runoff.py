import numpy as np
import os
import shutil
from multiprocessing import Pool
from multiprocessing import Process
from datetime import datetime, timedelta
import random
import netCDF4 as nc
import calendar

#external python codes
import params_case14 as pm
# ################################################
#
# make folders
# Prepare ensembles for DA
#
# ################################################
def read_nc(filepath):
        ds = nc.Dataset(filepath, mode='r')
        qr = ds.variables['Qr'][:]
        time = ds.variables['time']
        time_var = ds.variables['time'][:]
        time_units = time.units
        time_value = nc.num2date(time_var,time_units)
        start_time = time_value[0]
        end_time = time_value[-1]
        ds.close()
        qr_new = qr*86400 # unit convert
        return qr_new, start_time, end_time

def mv_ori_nc(new_path,bin_path,sim_path):
        os.makedirs(new_path, exist_ok=True)
        os.makedirs(bin_path, exist_ok=True)
        for filename in os.listdir(sim_path):
                if filename.endswith('Qr.nc'):
                        file_path = os.path.join(sim_path,filename)
                        mv_path = os.path.join(new_path,'Qr00.nc')
                        shutil.move(file_path,mv_path)
                        shutil.rmtree(sim_path)

def mv_ens_nc(new_path,bin_path,sim_path,ens_ind):
        os.makedirs(new_path, exist_ok=True)
        os.makedirs(bin_path, exist_ok=True)
        ens_str = '{:02d}'.format(ens_ind+1)
        for filename in os.listdir(sim_path):
                if filename.endswith('Qr.nc'):
                        file_path = os.path.join(sim_path,filename)
                        mv_path = os.path.join(new_path,'Qr' + ens_str + '.nc')
                        shutil.move(file_path,mv_path)
                        shutil.rmtree(sim_path)

def cal_error(ind, qr, error_use, new_path, mon_str):
        et = error_use[ind,:,:]
        et_use = np.broadcast_to(et[np.newaxis,:,:], qr.shape)
        result = et_use * qr
        roff32 = result.astype(np.float32)
        roff_flatten = roff32.flatten()
        file_size = os.path.getsize(new_path + 'Qr00.nc')
        if file_size > 30 * 5 * 24 * 1320 * 1500 * 4:
                filen = 'roff_mon' + mon_str + '_' + '{:02d}'.format(ind+1) + '.bin'
        else:
                filen = 'roff' + '{:02d}'.format(ind+1) + '.bin'
    
        file_name = new_path + filen
        with open(file_name, 'wb') as f:
                f.write(roff_flatten)
        return result


def pertubate_runoff(qr, mon, new_path):
        #runoff pertubation
        delta = 1-1/240
        et = np.full((pm.ens_mem(),np.shape(qr)[1],np.shape(qr)[2]),1.)
        a = np.random.normal(1,0.25,pm.ens_mem()*np.shape(qr)[1]*np.shape(qr)[2])
        std = np.sqrt((1-delta**2)*a**2)
        error = np.random.normal(1,std)
        error_use = error.reshape((pm.ens_mem(),np.shape(qr)[1],np.shape(qr)[2]))
        mon_str = '{:02d}'.format(mon+1)
        with Pool(pm.cpu_nums()) as p:
                new_runoff = p.starmap(cal_error, [(ind, qr, error_use, new_path, mon_str) for ind in range(pm.ens_mem())])
        return np.array(new_runoff)


def cal_runoff():
        ####### Main function to calculate runoff #############
        mode = pm.mode()
        force_mode = pm.force_mode(mode)
        bin_path = pm.runoff_dir() + '/tej/'
        mon = 0
        start_year,start_month,start_date,start_hour=pm.starttime()
        yyyy="%04d"%(start_year)
        if force_mode == 2:
                real_mode = pm.real_mode()
                if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                        # do DA for historical period
                        new_path = pm.ILS_dir() + '/his' + yyyy + '/'
                else:
                        # do DA for real time
                        new_path = pm.ILS_dir() + '/realtime/'
                #########################
                file_size = os.path.getsize(new_path + 'Qr00.nc')
                if file_size > 30*5*24*1320*1500*4:
                        # break it into smaller files
                        for mon in range(0,12):
                                ncpath= new_path + 'Qr00_mon' + '{:02d}'.format(mon+1) + '.nc'
                                qr, sta_time, end_time = read_nc(ncpath)
                                current_time = sta_time
                                time_ind = 0
                                result = pertubate_runoff(qr, mon, new_path)
                                while current_time <= end_time:
                                        std_filen = 'Runoff_' + current_time.strftime('%Y%m%d') + '000.bin'
                                        shutil.copy(pm.runoff_dir() + '/inp/' + std_filen, pm.runoff_dir() + '/tej/' + std_filen)
                                        for ens_ind in range(0,pm.ens_mem()):
                                                ens_str = '{:03d}'.format(ens_ind+1)
                                                if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                                                        qr_each = result[ens_ind,time_ind:time_ind+24,:,:]
                                                else:
                                                        qr_each = result[ens_ind,:,:]
                                                data32 = qr_each.astype(np.float32)
                                                data_flatten = data32.flatten()
                                                file_name = bin_path + 'Runoff_' + current_time.strftime('%Y%m%d') + ens_str + '.bin'
                                                with open(file_name, 'wb') as f:
                                                        f.write(data_flatten)
                                        time_ind = time_ind+24
                                        current_time += timedelta(days=1)

                ############### 
                else:
                        ncpath= new_path + 'Qr00.nc'
                        qr, sta_time, end_time = read_nc(ncpath)
                        new_runoff = pertubate_runoff(qr, mon, new_path)
                        time_ind = 0
                        current_time = sta_time
                        while current_time <= end_time:
                                for ens_ind in range(0,pm.ens_mem()):
                                        ens_str = '{:03d}'.format(ens_ind+1)
                                        if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                                                qr_each = new_runoff[ens_ind,time_ind:time_ind+24,:,:]
                                        else:
                                                qr_each = new_runoff[ens_ind,:,:,:]
                                        data32 = qr_each.astype(np.float32)
                                        data_flatten = data32.flatten()
                                        file_name = bin_path + 'Runoff_' + current_time.strftime('%Y%m%d') + '{:03d}'.format(ens_ind+1) + '.bin'
                                        with open(file_name, 'wb') as f:
                                                f.write(data_flatten)
                                time_ind = time_ind+24
                                current_time += timedelta(days=1)
                os.system('rm -rf '+new_path+'roff*')


def write_ori_bin(sim_path,new_path,bin_path,filen):
        ncpath= new_path + filen
        qr, sta_time, end_time = read_nc(ncpath)
        current_time = sta_time
        time_ind = 0
        real_mode = pm.real_mode()
        while current_time <= end_time:
                if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                        qr_each = qr[time_ind:time_ind+24,:,:]
                else:
                        qr_each = qr
                data32 = qr_each.astype(np.float32)
                data_flatten = data32.flatten()
                file_name = bin_path + 'Runoff_' + current_time.strftime('%Y%m%d') + '000.bin'
                with open(file_name, 'wb') as f:
                        f.write(data_flatten)
                time_ind = time_ind+24
                current_time += timedelta(days=1)


def write_ens_bin(sim_path,new_path,bin_path,filen,ens_ind):
        # ncpath= new_path + 'Qr' + '{:02d}'.format(ens_ind) + '.nc'
        ncpath= new_path + filen
        qr, sta_time, end_time = read_nc(ncpath)
        current_time = sta_time 
        time_ind = 0
        real_mode = pm.real_mode()
        while current_time <= end_time:
                if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                        qr_each = qr[time_ind:time_ind+24,:,:]
                else:
                        qr_each = qr
                data32 = qr_each.astype(np.float32)
                data_flatten = data32.flatten()
                file_name = bin_path + 'Runoff_' + current_time.strftime('%Y%m%d') + '{:03d}'.format(ens_ind+1) + '.bin'
                with open(file_name, 'wb') as f:
                        f.write(data_flatten)
                time_ind = time_ind+24
                current_time += timedelta(days=1)

def write_mon_bin(sim_path,bin_path,filen,ens_ind):
        ncpath= sim_path + filen
        start_time = datetime(*pm.starttime(DAhour=pm.DA_leading()))
        end_time   = datetime(*pm.endtime(DAhour=pm.DA_leading()))
        current_time = start_time
        time_ind = 0
        real_mode = pm.real_mode()
        while current_time <= end_time:
                yr_str  = current_time.strftime('%Y')
                day_str = current_time.strftime('%d')
                mon_str = current_time.strftime('%m')
                _,last_day = calendar.monthrange(int(yr_str),int(mon_str))
                if time_ind == 0:
                        qr, nc_statime, nc_endtime = read_nc(ncpath+mon_str+'.nc')
                if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                        qr_each = qr[time_ind:time_ind+24,:,:]
                else:
                        qr_each = qr
                data32 = qr_each.astype(np.float32)
                data_flatten = data32.flatten()
                file_name = bin_path + 'Runoff_' + current_time.strftime('%Y%m%d') + '{:03d}'.format(ens_ind+1) + '.bin'
                with open(file_name, 'wb') as f:
                        f.write(data_flatten)
                f.close()
                if current_time.day == last_day:
                        time_ind = 0
                else:
                        time_ind = time_ind+24
                current_time += timedelta(days=1)

def make_small_nc(sim_path, new_path, bin_path):
        # break it into smaller files
        file_size = os.path.getsize(new_path+'Qr00.nc')
        if file_size > 30*5*24*1320*1500*4:
                shutil.copy(new_path+'Qr00.nc', new_path+'Qr.nc')
                os.system('cdo splitmon '+new_path+'Qr.nc Qr00_mon')
                for mon in range(0,12):
                        mon_str = '{:02d}'.format(mon+1)
                        write_ori_bin(sim_path,new_path,bin_path,'Qr00_mon'+mon_str+'.nc')
        else:
                write_ori_bin(sim_path,new_path,bin_path,'Qr00.nc')
        
def process_ensemble(ens_ind, new_path, bin_path):
        ens_str = '{:02d}'.format(ens_ind+1)
        sim_path = pm.ILS_dir() + '/wlv' + ens_str + '/'
        if not os.path.exists(sim_path):
                print('ILS ensembles failed...')
        else:
                mv_ens_nc(new_path,bin_path,sim_path,ens_ind)
        file_size = os.path.getsize(new_path+'Qr'+ens_str+'.nc')
        if file_size > 30*5*24*1320*1500*4:
                shutil.copy(new_path+'Qr'+ens_str+'.nc', new_path+'Qr_mon'+ens_str+'.nc')
                os.system('cdo splitmon '+new_path+'Qr_mon'+ens_str+'.nc Qr'+ens_str+'_mon')
                for mon in range(0,12):
                        mon_str = '{:02d}'.format(mon+1)
                write_ens_bin(sim_path,new_path,bin_path,'Qr'+ens_str+'_mon'+mon_str+'.nc',ens_ind)
        else:   
                write_ens_bin(sim_path,new_path,bin_path,'Qr' + '{:02d}'.format(ens_ind+1) + '.nc',ens_ind)

def process_mon_ensemble(ens_ind, file_path, bin_path):
        ens_str = '{:02d}'.format(ens_ind+1)
        sim_path = pm.ILS_dir() + '/his2024/'
        write_mon_bin(sim_path,bin_path,'Qr'+'{:02d}'.format(ens_ind+1)+'_mon',ens_ind)

def wrapper(ens_ind, new_path, bin_path):
        # process_ensemble(ens_ind, new_path, bin_path)
        process_mon_ensemble(ens_ind, new_path, bin_path)

def ensembles():
        start_year,start_month,start_date,start_hour=pm.starttime()
        yyyy="%04d"%(start_year)
        mode = pm.mode()
        force_mode = pm.force_mode(mode)
        ens_num = pm.ens_mem()
        bin_path = pm.runoff_dir() + '/ils/'
        real_mode = pm.real_mode()
        ##############################
        # do DA for historical period
        if (real_mode == 0) | (real_mode == 1) | (real_mode == 2):
                new_path = pm.ILS_dir() + '/his' + yyyy + '/'
                os.makedirs(new_path, exist_ok=True)
                sim_path = pm.ILS_dir() + '/sim'+yyyy+'/'
                # non-pertubation data
                if not os.path.exists(sim_path):
                        print('ILS ensembles failed...')
                else:
                        mv_ori_nc(new_path,bin_path,sim_path)
                ### already split
                # make_small_nc(sim_path, new_path, bin_path) # generate cama/inp file
                if force_mode == 1:  
                        pool = Pool(processes=pm.cpu_nums())    
                        pool.starmap(wrapper, [(ens_ind, new_path, bin_path) for ens_ind in range(0,ens_num)])
                        pool.close()
                        pool.join()
        ###############################
        # do DA for real time
        else:
                new_path = pm.ILS_dir() + '/realtime/'
                os.makedirs(new_path, exist_ok=True)
                sim_path = pm.ILS_dir() + '/real/'
                # runoff pertubation
                if not os.path.exists(sim_path):
                        print('ILS ensembles failed...')
                else:
                        mv_ori_nc(new_path,bin_path,sim_path)
                make_small_nc(sim_path, new_path, bin_path) # generate cama/inp file
                if force_mode == 1:  
                        # rainfall pertubation
                        pool = Pool(processes=pm.cpu_nums())
                        pool.starmap(wrapper, [(ens_ind, new_path, bin_path) for ens_ind in range(0,ens_num)])
                        pool.close()
                        pool.join()
        #####################################


