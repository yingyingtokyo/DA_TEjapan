
def draw_lines(varname,lat_valid_plot, lon_valid_plot, assim_num_plot,obs_station,da_station,sim_station,savename,new_row,alabel,syyyy,tstart,time_len,simmean,pm_name,rainf,lon,lat,pm_path,expname,dahour_str,exp_plotdir):
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
    valid_indices = [i for i in range(len(lon_valid_plot)) if not np.all(np.isnan(obs_station[tstart:, i]))]
    if "Outflow" in varname:
        # to remove the very small river stations
        valid_indices = [i for i in range(len(lon_valid_plot)) if not (np.all(np.isnan(obs_station[tstart:, i])) or np.all(obs_station[tstart:, i] < 5))]
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
        station_flood_time = flood_time_valid[i,:]
        if "WSE" in varname:
            sim_plot = sim_station[tstart:,:,i]-mean_station[i]
            da_plot  = da_station[tstart:,:,i]-mean_station[i]
        else:
            sim_plot = sim_station[tstart:,:,i]
            da_plot  = da_station[tstart:,:,i]
        # RMSE
        rmse_da = cal_rmse(obs_plot,np.nanmean(da_plot,axis=1),varname,time_len)
        kge_da = cal_kge(obs_plot,np.nanmean(da_plot,axis=1),varname)
        nse_da = cal_nse(obs_plot,np.nanmean(da_plot,axis=1),varname)
        rmse_da_str = '{:.2f}'.format(rmse_da)
        # rmse_sim= cal_rmse(obs_plot,np.nanmean(sim_plot,axis=1),varname,time_len)  # compare with the average
        rmse_sim= cal_rmse(obs_plot,sim_plot[:,-1],varname,time_len)  # compare with open loop
        kge_sim= cal_kge(obs_plot,sim_plot[:,-1],varname,time_len)  # compare with open loop
        nse_sim= cal_nse(obs_plot,sim_plot[:,-1],varname,time_len)  # compare with open loop
        rmse_sim_str= '{:.2f}'.format(rmse_sim)
        nse_point = nse_point + [[nse_da,nse_sim]]
        kge_point = kge_point + [[kge_da,kge_sim]]
        if station_flood_time[0]<-900:
            # no flood event
            improve_rate[ind] = np.nan
            continue
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
        ax2.bar(time[tstart:],prep_station[tstart:,i],color='gray',alpha=0.8)
        ax2.set_ylim(np.nanmax(prep_station[tstart:,i])+5,0)
        ax2.spines['right'].set_color('gray')
        ax2.yaxis.label.set_color('gray')
        ax2.tick_params(axis='y',colors='gray') 
        # floodtime
        # axes[row,col].axvline(x=station_flood_time[0],color='gray',linestyle='--',lw=1.)
        # axes[row,col].axvline(x=station_flood_time[1],color='gray',linestyle='--',lw=1.)
        # simulations
        sim_point_ind = sim_point_ind + [sim_plot]
        axes[row,col].plot(time_plot,sim_plot,'b-',linewidth=0.2,alpha=0.8)
        # online DA
        da_point_ind  = da_point_ind + [da_plot]
        axes[row,col].plot(time_plot,da_plot,'r-',linewidth=1.25)
        # print(np.unique(da_plot))
        axes[row,col].set_xticks(np.arange(tstart+15,time_len,24))
        axes[row,col].set_xlim(tstart,time_len)
        # open loop
        axes[row,col].plot(time_plot,sim_plot[:,-1],'b-',linewidth=1.5)
        # axes[row,col].plot(time_plot,np.nanmean(sim_plot,axis=1),'b-',linewidth=1.5)
        # obs
        obs_point_ind = obs_point_ind + [obs_plot]
        axes[row,col].plot(time_plot,obs_plot,color='green',linewidth=2.5,linestyle='--')
        lat_str = '{:.2f}'.format(lat[int(lat_valid_plot[i])-1])
        # lat_str = '{:.2f}'.format(int(lat_valid_plot[i])-1)
        lon_str = '{:.2f}'.format(lon[int(lon_valid_plot[i])-1])
        # lon_str = '{:.2f}'.format(int(lon_valid_plot[i])-1)
        if rmse_da < rmse_sim:
            improve_rate[ind] = 1
        rmse_all = rmse_all+ [plot_idx,rmse_da,rmse_sim,kge_da,kge_sim,int(lat_valid_plot[i])-1,int(lon_valid_plot[i])-1,int(assim_num_plot[i])]
        axes[row,col].set_title(alabel[plot_idx],loc='left',fontsize=12)
        axes[row,col].text(0.02,0.95,'('+lat_str+','+lon_str+')\n' + 'obs number: '+str(int(assim_num_plot[i])),transform=axes[row,col].transAxes,fontsize=12, c='black',fontweight='bold', va='top', ha='left')
        axes[row,col].text(0.9,0.95,'RMSE: '+rmse_da_str,transform=axes[row,col].transAxes,fontsize=12, c='red',fontweight='bold', va='top', ha='right')
        axes[row,col].text(0.9,0.85,'RMSE: '+rmse_sim_str,transform=axes[row,col].transAxes,fontsize=12, c='blue',fontweight='bold', va='top', ha='right')
        if col == 0:
            axes[row,col].set_ylabel(varname)
        if col == cols-1:
            ax2.set_ylabel('Rain (mm/h)')
        if 'ils01' in pm_path or 'round' or 'msm' or 'all' in pm_path:
            if 'WSE' in varname:
                if row == 7: # parmas01
                    dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='24h')
                    axes[row,col].set_xticklabels(dates,rotation=25)
                else:
                    axes[row,col].set_xticklabels([],rotation=15)
            else:
                    if row == 5:
                        dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='24h')
                        axes[row,col].set_xticklabels(dates,rotation=25)
                    else:
                        axes[row,col].set_xticklabels([],rotation=15)     
        if 'ils03' in pm_path:
            if 'WSE' in varname:
                if row == rows-1: # parmas03
                    dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='24h')
                    axes[row,col].set_xticklabels(dates,rotation=25)
                else:
                    axes[row,col].set_xticklabels([],rotation=15)
            else:
                    if row == 5:
                        dates = pd.date_range(start='2019-10-11 00:00',end='2019-10-17 08:00',freq='24h')
                        axes[row,col].set_xticklabels(dates,rotation=25)
                    else:
                        axes[row,col].set_xticklabels([],rotation=15) 
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
    nse_point = np.array(nse_point)
    kge_point = np.array(kge_point)
    nse_da   = np.shape(np.where(nse_point[:,0]>0.5))[1]
    nse_da_good  = nse_da/np.shape(np.where(~np.isnan(improve_rate)))[1]*100
    nse_sim  = np.shape(np.where(nse_point[:,1]>0.5))[1]
    nse_sim_good = nse_sim/np.shape(np.where(~np.isnan(improve_rate)))[1]*100
    print(np.nansum(improve_rate),np.shape(np.where(~np.isnan(improve_rate)))[1])
    print('How many stations improve: ','{:.2f}'.format(improve),'%')
    print('NSE_da_good: ','{:.2f}'.format(nse_da_good),'%   ','NSE_sim_good: ','{:.2f}'.format(nse_sim_good),'%')
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
            plt.savefig(exp_plotdir+'2.cross_valcase'+savename+pm_path[-4]+'.jpg', format='jpg',dpi=600)
        else:
            np.save(pm.plot_dir()+'/exp_'+expname+dahour_str+'/rmse'+savename+'.npy',np.array(rmse_all))
            plt.savefig(exp_plotdir+'2.cross_validation'+savename+'.jpg', format='jpg',dpi=600)
    plt.show()
    plt.close()
    return np.array(lon_point), np.array(lat_point), np.array(lon_point_ind), np.array(lat_point_ind), np.array(obs_point_ind), np.array(sim_point_ind), np.array(da_point_ind), np.array(num_point), new_row





def draw_each_figure(obs_lat,obs_lon,varname):
    situ_rmse  = np.full(np.shape(obs_lat)[0],np.nan)
    situ_nrmse = np.full(np.shape(obs_lat)[0],np.nan)
    situ_mse = np.full((np.shape(obs_lat)[0],2),np.nan)
    situ_kge   = np.full((np.shape(obs_lat)[0],2),np.nan)
    situ_nse   = np.full((np.shape(obs_lat)[0],2),np.nan)
    lat_plot   = np.full(np.shape(obs_lat)[0],1)
    lon_plot   = np.full(np.shape(obs_lat)[0],1)
    assim_plot = np.full(np.shape(obs_lat)[0],1)
    situ_imp   = np.full(np.shape(obs_lat)[0],np.nan)
    loc_row = np.full(np.shape(obs_lat)[0],np.nan)
    if varname == 'rivdph':
        typename = 'rivdph'
    else:
        typename = 'outflw'
    if 'all' in pm_path:
        if 'case' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+pm_name[-5:]+'/'+typename+'/all/'
        else:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+typename+'/'
    else:
        if 'round' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/round'+pm_path[-4]+'/'+typename+'/'
        elif 'case' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+pm_name[-5:]+'/'+typename+'/valid/'
        elif 'ens' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/MEPS/'+typename+'/'
        else:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/'+typename+'/'
    # DA
    filemean = filepath+'meanA.bin'   # compare to the average of ensemble
    da_mean = read_bin_mean(filemean)  # compare with the open loop result
    # sim
    filemean = filepath+'meanC.bin'  # compare with average of ensembles
    # filemean = filepath+'oriC.bin' # compare with original simulations
    sim_mean = read_bin_mean(filemean)  # compare with the open loop result
    # if '00' in filename:
        # sim_mean = read_nc(filename,4)
    if varname == 'HQ':
        if 'round' in pm_path:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/round'+pm_path[-4]+'/rivdph/'
        else:
            filepath = inputdir+'/exp_'+expname+dahour_str+'/rivdph/'
        da_dis_mean = read_bin_mean(filepath+'meanA.bin')
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
        if 'case1' not in pm_name:
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
            if np.isnan(mse_sim) | np.isnan(mse_da) | np.all(np.isnan(loc_da)) | (np.nanmean(loc_obs)<obs_range): #| (np.nanmean(loc_sim)<obs_range) | (np.nanmean(loc_da)<obs_range):
                continue
            rmse = np.sqrt(mse_sim) - np.sqrt(mse_da)
            if (np.shape(loc_obs)[0]<25) | (np.shape(loc_sim)[0]<25) | (np.shape(loc_da)[0]<25): #| (np.abs(elvmean[loc_lat,loc_lon]-simmean[loc_lat,loc_lon])>5):
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
                    _,_,situ_kge[loc,0],situ_nse[loc,0],_  = cal_single_metric(loc_obs,loc_da)
                    _,_,situ_kge[loc,1],situ_nse[loc,1],_  = cal_single_metric(loc_obs,loc_sim)
                imp = rmse/np.sqrt(mse_sim)            
                # exclude lake region
                if rmse<-0.:
                    if np.abs(elvmean[loc_lat,loc_lon]-simmean[loc_lat,loc_lon])>3:
                        rmse = np.nan
                        nrmse = np.nan
                        imp = np.nan 
                        situ_kge[loc,:] = np.nan
                        situ_nse[loc,:] = np.nan
                    # else:
                    #     print(loc,loc_lat,loc_lon,np.nanmean(loc_obs))
            # store nrmse, lat, loc
            situ_rmse[loc]=rmse
            situ_mse[loc,0]=np.sqrt(mse_da)
            situ_mse[loc,1]=np.sqrt(mse_sim)
            situ_nrmse[loc]=nrmse
            situ_imp[loc]=imp
            lat_plot[loc]=int(loc_lat)
            lon_plot[loc]=int(loc_lon)
            assim_plot[loc] = int(loc_assim)
            loc_row[loc]=loc
    return situ_rmse,situ_mse,situ_nrmse,situ_imp,loc_row,lat_plot,lon_plot,assim_plot,situ_kge,situ_nse
