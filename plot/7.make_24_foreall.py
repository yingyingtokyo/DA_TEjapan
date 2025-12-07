import numpy as np
from datetime import datetime, timedelta
import netCDF4 as nc
import sys
from multiprocessing import Pool
params_path = '/data42/yingying/HydroDA/src'
#external python codes
sys.path.append(params_path)
import params_ils03 as pm

camadir  = pm.CaMa_dir()
inputdir = pm.plot_dir()
expname  = pm.runname(pm.mode())
dahour   = pm.dahour()
dahour_str = '{:02d}'.format(dahour)
exp_plotdir = inputdir + '/exp_' + expname + dahour_str + '/'
obsdir = pm.obs_dir(pm.real_mode())

# simulation start from 20191007 but analysis start from 20191010
start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
syyyy='%04d' % (start_year)
if syyyy == '2019':
    tstart = int(24*3/dahour)
else:
    tstart = 0
smm='%02d' % (start_month)
sdd='%02d' % (start_date)
shh='%02d' % (start_hour)
stime = datetime(start_year,start_month,start_date,start_hour)
end_year,end_month,end_date,end_hour=pm.endtime() # Start year month date
eyyyy='%04d' % (end_year)
emm='%02d' % (end_month)
edd='%02d' % (end_date)
ehh='%02d' % (end_hour)
etime = datetime(end_year,end_month,end_date,end_hour)
time_step = timedelta(hours=dahour)
time_len = int((etime-stime).total_seconds() / 3600 / dahour)
lon_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lon_point_ind.npy')
lat_point_ind = np.load(pm.plot_dir()+'/exp_'+expname+dahour_str+'/lat_point_ind.npy')

def read_nc(fname,layer):
    nf = nc.Dataset(fname,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    lat = np.array(nf.variables[varname[1]][:])
    lon = np.array(nf.variables[varname[2]][:])
    var = np.array(nf.variables[varname[layer]][:])
    var = np.array(np.where(var>10**7,np.nan,var))
    return var   

# read 24h data
def read_24h_nc(enum,varname):
    ind = 0
    current_time = stime
    filepath = pm.outdir()+'/CaMa_out_'+expname+dahour_str+'/'
    end_date = etime-timedelta(days=1)
    ens = '{:03d}'.format(enum+1)
    wlv_da_each = np.full((time_len-int(24/dahour),np.shape(lon_point_ind)[0],24),np.nan)  # 24 hours forecast
    wlv_sim_each = np.full((time_len-int(24/dahour),np.shape(lon_point_ind)[0],24),np.nan)  # 24 hours forecast
    while current_time < end_date:
        ctime = current_time.strftime('%Y%m%d%H')
        args = [(enum,varname) for enum in range(pm.ens_mem())]
        nyyyy = current_time.strftime('%Y')
        data_wlv_da  = np.full((np.shape(lon_point_ind)[0],24),np.nan)
        data_wlv_sim = np.full((np.shape(lon_point_ind)[0],24),np.nan)
        fname_da = filepath+ctime+'A'+ens+'/'+'o_'+varname+nyyyy+'.nc'
        fname_sim = filepath+ctime+'C'+ens+'/'+'o_'+varname+nyyyy+'.nc'
        data_file_da  = read_nc(fname_da,3)
        data_file_sim = read_nc(fname_sim,3)
        for i in range(0,np.shape(lon_point_ind)[0]):
            data_wlv_da[i,:]   = data_file_da[:,lat_point_ind[i],lon_point_ind[i]]
            data_wlv_sim[i,:]  = data_file_sim[:,lat_point_ind[i],lon_point_ind[i]]
        wlv_da_each[ind,:,:]  = data_wlv_da
        wlv_sim_each[ind,:,:] = data_wlv_sim
        current_time += time_step
        ind = ind + 1
    return wlv_da_each,wlv_sim_each

def parallel(varname):
    args = [(enum,varname) for enum in range(pm.ens_mem())]
    with Pool(processes=pm.cpu_nums()) as p:
        results = p.starmap(read_24h_nc,args)
    wlv_da = np.full((pm.ens_mem(),time_len-int(24/dahour),np.shape(lon_point_ind)[0],24),np.nan)  # 24 hours forecast
    wlv_sim = np.full((pm.ens_mem(),time_len-int(24/dahour),np.shape(lon_point_ind)[0],24),np.nan)  # 24 hours forecast
    for i,(da,sim) in enumerate(results):
        wlv_da[i]=da
        wlv_sim[i]=sim
    savepath = pm.plot_dir()+'/exp_'+expname+dahour_str+'/'+varname+'/'
    np.save(savepath+'da_24fore.npy',wlv_da)
    np.save(savepath+'sim_24fore.npy',wlv_sim)
parallel('rivdph')
parallel('outflw')
