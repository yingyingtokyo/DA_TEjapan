#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#libralies
import os
import itertools
import numpy as np
import sys
import errno
from multiprocessing import Pool
from multiprocessing import Process
import datetime
import functools
import numpy.random as rd
import os.path
import datetime as dt
import glob
import shutil
from numpy import ma
import random
import re
import calendar
import math
import netCDF4 as nc

#external python codes
import params_round5 as pm
########################
## main control function
############
def main_action():
    """Main setting for experiment"""
    mode=pm.mode()
    run_name=pm.runname(pm.mode)
    # name of restart file after spin-up
    mm_restart = "%02d"%(5)
    dd_restart = "%02d"%(1)
    # DA setting
    # experiment time
    start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
    end_year,end_month,end_date,end_hour=pm.endtime() # End year month date
    start_dt=datetime.datetime(start_year,start_month,start_date,start_hour)
    end_dt=datetime.datetime(end_year,end_month,end_date,end_hour)
    # how many hours
    hour_count=int((end_dt-start_dt).total_seconds()/3600/pm.dahour())


    """Main function to prepare experiment settings"""
    print('Exp:',pm.runname(pm.mode()))
    initial_setting(mm_restart,dd_restart)

    # spin-up
    #print ("spin up simulation")
    #spin_up()
    #exit()

    # calculate average of simulation for DA
    calc_mean()


    for hour in range(0,hour_count): 
        running_dt=start_dt+datetime.timedelta(hours=pm.dahour()*hour)
        yyyy='%04d' % (running_dt.year)
        mm='%02d' % (running_dt.month)
        dd='%02d' % (running_dt.day)
        hh='%02d' % (running_dt.hour)
        # CaMa simulation
        one_hour_loop(yyyy,mm,dd,hh)

    # clean all intermediate files
    if pm.output_er()==1:
        os.system("rm -Rf "+pm.expout()+"/"+pm.CaMa_out_dir()+"/"+yyyy+"*")

def make_initial_infl():
    parm_infl=np.ones((1320,1500),dtype=np.float32)
    start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
    yyyy='%04d' % (start_year)
    mm='%02d' % (start_month)
    dd='%02d' % (start_date)
    hh='%02d' % (start_hour)
    file_path = pm.outdir()+pm.CaMa_out_dir()+"/inflation"
    os.makedirs(file_path,exist_ok=True)
    parm_infl.tofile(file_path+"/parm_infl"+yyyy+mm+dd+hh+".bin")

############
def initial_setting(mm_restart,dd_restart):
    start_year,start_month,start_date,start_hour=pm.starttime()
    yyyy="%04d"%(start_year)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    hh="%02d"%(start_hour)
    inputlist=[] # parallel
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        spinup_open="%04d%2d%02dC%03d"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date(),num)
        os.makedirs(pm.expout()+"/"+pm.CaMa_out_dir()+'/'+spinup_open,exist_ok=True)
        # Don't change this part, the filename here is fixed
        os.system("cp "+pm.expout()+'/'+pm.CaMa_out()+"/"+spinup_open+"/restart2019"+mm_restart+dd_restart+".bin "+pm.expout()+'/'+pm.CaMa_out_dir()+"/"+spinup_open+"/restart"+yyyy+mm+dd+".bin")
        os.system("cp "+pm.expout()+'/'+pm.CaMa_out()+"/"+spinup_open+"/restart2019"+mm_restart+dd_restart+".bin "+pm.expout()+'/'+pm.CaMa_in_dir()+"/restart/open/restart"+yyyy+mm+dd+hh+"C"+numch+".bin")
        os.system("cp "+pm.expout()+'/'+pm.CaMa_out()+"/"+spinup_open+"/restart2019"+mm_restart+dd_restart+".bin "+pm.expout()+'/'+pm.CaMa_in_dir()+"/restart/assim/restart"+yyyy+mm+dd+hh+"A"+numch+".bin")
    return 0


############
## main program functions
############
"""spin up code"""
def spinup_loop(inputlist):
    # Run spinup simulation
    run_name=pm.runname(pm.mode)
    yyyy, loop, ens_num = inputlist
    dir2=pm.CaMa_dir()
    cpunums=pm.cpu_nums()
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    mapname=pm.mapname()
    cal=pm.calibrate()
    print("%s for %03d"%(loop,int(ens_num)))
    os.system("source "+str(pm.DA_dir())+"/src/spin_up.sh "+str(yyyy)+" "+str(loop)+" "+ens_num+" "+dir2+" "+str(cpunums)+" "+str(run_name)+" "+str(exp_dir)+" "+str(mapname)+" "+str(cal)+" "+pm.CaMa_out_dir()+" "+pm.CaMa_in_dir())
    return 0

def spin_up(): #used
    # run spin up simulation
    # 1 year spin up for calculating initial value
    # one simulation for true
    # ensmble simulation for open
    cpunums = pm.cpu_nums()
    inputlist=[]
    yyyy=str(pm.spinup_end_year())
    for ens_num in np.arange(1,pm.ens_mem()+1):
        inputlist.append([yyyy,"open",'%03d'%ens_num])
    # Run spinup simulations
    p=Pool(pm.para_nums())
    p.map(spinup_loop,inputlist)
    p.terminate()
    print("======================= end spinup ==========================")
    return 0

###########################
"""CaMa simulation code"""

def one_hour_sim(inputlist):
    yyyy,mm,dd,hh,ens_num,looptype = inputlist
    da_ehour_str = '{:02d}'.format(pm.da_ehour()) #DA frequency
    run_name=pm.runname(pm.mode())
    # program for running one day model
    bef_dt=datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))-datetime.timedelta(hours=pm.dahour())
    bef_yyyy='%04d' %bef_dt.year
    bef_mm='%02d' %bef_dt.month
    bef_dd='%02d' %bef_dt.day
    bef_hh='%02d' %bef_dt.hour

    print ("oneday loop for",yyyy,mm,dd,hh,ens_num,looptype)
    dir2=pm.CaMa_dir()
    lead_str = pm.DA_leading() #DA leading time

    print (yyyy+" "+mm+" "+dd+" "+hh+" "+ens_num+" "+dir2+" "+looptype)
    cpunums = pm.cpu_nums() # Moiz: OMP takes integers
    exp_dir=pm.expout() #pm.DA_dir()+"/out/"+pm.experiment()
    mapname=pm.mapname()
    cal=pm.calibrate()
    DA_dir=pm.DA_dir()
    os.system("source "+pm.DA_dir()+"/src/onehour_sim.sh "+yyyy+" "+mm+" "+dd+" "+hh+" "+ens_num+" "+dir2+" "+looptype+" "+str(cpunums)+" "+str(run_name)+" "+str(exp_dir)+" "+str(mapname)+" "+str(cal)+" "+DA_dir+" "+pm.input(pm.mode())+" "+da_ehour_str+' '+pm.CaMa_out_dir()+" "+pm.CaMa_in_dir()+" "+lead_str+" "+str(pm.val_mode())) #+" "+str(DT)) #+str(corrupt)+" "
    return 0

def one_hour_loop(yyyy,mm,dd,hh):
    print ("===== start loop of "+yyyy+" "+mm+" "+dd+" "+hh+" =====")
    # set up the initial influtation
    make_initial_infl()

    # Corrupted Simulation (Open Loop) #####
    if pm.run_flag() == 0:
        ODM_inputlist=[]
        # set for ensemble simulations
        ens_num=1
        for ens_num in np.arange(1,pm.ens_mem()+1):
            ODM_inputlist.append([yyyy,mm,dd,hh,'%03d'%ens_num,"open"])
        # Run CaMa-Flood Model (ensemble simulations),read restart from ori result
        p=Pool(pm.para_nums())
        p.map(one_hour_sim,ODM_inputlist)
        p.terminate()
        #copy ensemble corrupted simulation forecasts from CaMa_out (xa)
        cprestart=[]
        for ens_num in np.arange(1,pm.ens_mem()+1):
            cprestart.append([yyyy,mm,dd,hh,ens_num])
        p=Pool(pm.para_nums())
        p.map(copy_corrupted_sfcelv,cprestart)
        p.terminate()

        p=Pool(pm.para_nums())
        p.map(copy_corrupted_restart,cprestart)
        p.terminate()

    # Assimilated Simulation ###################
    ODM_inputlist=[]
    # set for ensemble simulations, read restart from DA result
    for ens_num in np.arange(1,pm.ens_mem()+1):
        ODM_inputlist.append([yyyy,mm,dd,hh,'%03d'%ens_num,"assim"])
    # Run CaMa-Flood Model (ensemble simulations)
    p=Pool(pm.para_nums())
    p.map(one_hour_sim,ODM_inputlist)
    p.terminate()

    # make forecasted value for assimilated simulation
    # do assimilation (LETKF)
    data_assim(yyyy,mm,dd,hh)
    mkrestart=[]
    for ens_num in np.arange(1,pm.ens_mem()+1):
        mkrestart.append([yyyy,mm,dd,hh,"assim",'%03d'%ens_num])
    # Modify the restart file
    p=Pool(pm.para_nums())
    p.map(make_restart,mkrestart)
    p.terminate()
    # store river variable files
    store_out(yyyy,mm,dd,hh)
    # clean files
    if pm.output_er()==1:
        bef_dt=datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))-datetime.timedelta(hours=pm.dahour())
        bef_yyyy='%04d' %bef_dt.year
        bef_mm='%02d' %bef_dt.month
        bef_dd='%02d' %bef_dt.day
        bef_hh='%02d' %bef_dt.hour
        if bef_dt != datetime.datetime(bef_dt.year,12,31,23): # to keep restart file in last hour of the year
            if pm.val_mode()==0:
                os.system("rm -Rf "+pm.expout()+'/'+pm.CaMa_out_dir()+"/"+bef_yyyy+bef_mm+bef_dd+bef_hh+"*")

###########################
def read_nc(file):
    # 20190301-20190501
    # average of 201904
    # check if the existence of data
    if os.path.exists(file):
        nf = nc.Dataset(file,'r')
        varname = nf.variables.keys()
        varname = list(varname)
        print(varname)
        #[time,lat,lon]
        var = np.array(nf.variables[varname[4]])
        var = np.where(var>10000,np.nan,var)
        var = np.where(var<-100,np.nan,var)
        print(var[0,426,1016])
        var_data=var[0,:,:]
    return var_data
########################### # modified to run paralle @Menaka 

"""copy files for restarting"""
def copy_corrupted_sfcelv(inputlist):
    dahour_str = '{:02d}'.format(pm.dahour()) #DA frequency
    yyyy,mm,dd,hh,num = inputlist
    numch='%03d'%num
    fname=pm.expout()+'/'+pm.CaMa_out_dir()+"/"+yyyy+mm+dd+hh+"C"+numch+"/sfcelv"+yyyy+".bin"
    os.system("cp "+fname+" "+pm.expout()+"/assim_"+(pm.runname(pm.mode()))+dahour_str+"/ens_xa/open/"+yyyy+mm+dd+hh+"_"+numch+"_xa.bin")
    return 0

def copy_corrupted_restart(inputlist):
    yyyy,mm,dd,hh,num = inputlist
    nxt_day = datetime.datetime(int(yyyy),int(mm),int(dd),int(hh)) + datetime.timedelta(hours=pm.dahour())
    n_yyyy='%04d' % (nxt_day.year)
    n_mm='%02d' % (nxt_day.month)
    n_dd='%02d' % (nxt_day.day)
    n_hh='%02d' % (nxt_day.hour)
    numch='%03d'%num
    fname=pm.expout()+'/'+pm.CaMa_out_dir()+"/"+yyyy+mm+dd+hh+"C"+numch+"/restart"+n_yyyy+n_mm+n_dd+n_hh+".bin"
    print('copy',fname)
    print ("copy restart",n_yyyy,n_mm,n_dd,n_hh,"C"+numch)
    os.system("cp "+fname+" "+pm.expout()+'/'+pm.CaMa_in_dir()+"/restart/open/restart"+n_yyyy+n_mm+n_dd+n_hh+"C"+numch+".bin") ## CaMa-Flood v4.1
    return 0

###########################
def data_assim(yyyy,mm,dd,hh): # new data assimilation function (2020/05/18)
    dahour_str = '{:02d}'.format(pm.dahour()) #DA frequency
    parallels="%d"%(pm.para_nums()*pm.cpu_nums())
    os.environ['OMP_NUM_THREADS']=parallels
    exp_dir=pm.expout()
    thisday=datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))
    nxt_day=thisday+datetime.timedelta(hours=pm.dahour())
    nyear=nxt_day.year
    nmon=nxt_day.month
    nday=nxt_day.day
    nhour=nxt_day.hour
    print ('****DA****')
    nxtyyyymmddhh="%04d%02d%02d%02d"%(nyear,nmon,nday,nhour)
    os.system(pm.DA_dir()+"/src/data_assim "+yyyy+mm+dd+hh+" "+str(pm.patch_size())+" "+str(pm.ens_mem(pm.mode()))+" "+nxtyyyymmddhh+" "+pm.CaMa_dir()+" "+str(pm.thersold())+" "+exp_dir+" "+pm.DA_dir()+" "+pm.patch_filedir()+" "+" "+pm.obs_dir()+" "+pm.alloc_dir()+" "+str(pm.rho())+" "+str(pm.sigma_b())+" "+str(pm.conflag())+" "+(pm.calibrate())+" "+(pm.runname(pm.mode()))+" "+pm.CaMa_out_dir()+" "+pm.CaMa_in_dir()+" assim"+" "+str(pm.dahour())+" "+dahour_str+" "+str(pm.val_mode())+" "+str(pm.val_round(pm.val_mode())))  # online DA
    os.system(pm.DA_dir()+"/src/data_assim "+yyyy+mm+dd+hh+" "+str(pm.patch_size())+" "+str(pm.ens_mem(pm.mode()))+" "+nxtyyyymmddhh+" "+pm.CaMa_dir()+" "+str(pm.thersold())+" "+exp_dir+" "+pm.DA_dir()+" "+pm.patch_filedir()+" "+" "+pm.obs_dir()+" "+pm.alloc_dir()+" "+str(pm.rho())+" "+str(pm.sigma_b())+" "+str(pm.conflag())+" "+(pm.calibrate())+" "+(pm.runname(pm.mode()))+" "+pm.CaMa_out_dir()+" "+pm.CaMa_in_dir()+" open"+" "+str(pm.dahour())+" "+dahour_str+" "+str(pm.val_mode())+' '+str(pm.val_round(pm.val_mode())))  # offline DA
    return 0
###########################
def store_out(yyyy,mm,dd,hh): # update on 2023/05/30
    # program for storing data #
    dahour_str = '{:02d}'.format(pm.dahour()) #DA frequency
    listCA = ["open","assim"]
    #==========================#
    for looptype in listCA: # update on 2023/05/30
        if looptype == "open":
            CA = "C"
        else:
            CA = "A"

        for num in np.arange(1,pm.ens_mem()+1):
            numch = '%03d' % num 
            for var in ["rivdph","sfcelv","outflw","storge"]: 
                shutil.copy(pm.expout()+'/'+pm.CaMa_out_dir()+"/"+yyyy+mm+dd+hh+CA+numch+"/"+var+yyyy+".bin",pm.expout()+"/assim_"+(pm.runname(pm.mode()))+dahour_str+"/"+var+"/"+looptype+"/"+var+yyyy+mm+dd+hh+"_"+numch+".bin")

###########################
def make_restart(inputlist):
    dahour_str = '{:02d}'.format(pm.dahour()) #DA frequency
    # Unpack input parameters
    yyyy, mm, dd, hh, loop, numch  = inputlist
    print ("finish assimilating")

    # get the date of one day before
    nowdate=dt.datetime(int(yyyy),int(mm),int(dd),int(hh))
    befdate=nowdate-dt.timedelta(hours=pm.dahour())
    yyyy_b="%04d"%befdate.year
    mm_b="%02d"%befdate.month
    dd_b="%02d"%befdate.day
    hh_b="%02d"%befdate.hour

    # get the date of one day after
    nextdate=nowdate+dt.timedelta(hours=pm.dahour())
    yyyy_n="%04d"%nextdate.year
    mm_n="%02d"%nextdate.month
    dd_n="%02d"%nextdate.day
    hh_n="%02d"%nextdate.hour

    # calculate other variables from water storage
    exp_dir=pm.expout()
    os.system(pm.DA_dir()+"/src/make_restart "+yyyy+mm+dd+hh+" "+yyyy_b+mm_b+dd_b+hh_b+" "+yyyy_n+mm_n+dd_n+hh_n+" "+loop+" "+pm.CaMa_dir()+" "+pm.mapname()+" "+numch+" "+exp_dir+" "+pm.CaMa_out_dir()+" "+pm.CaMa_in_dir()+" "+(pm.runname(pm.mode()))+" "+dahour_str)

    print ("finish restarting",numch)
###########################
def read_sim_data(ens):
    year=pm.spinup_end_year()
    month=pm.spinup_end_month()
    date=pm.spinup_end_date()
    yyyy="%04d"%(year)
    mm="%02d"%(month)
    dd="%02d"%(date)
    ens_str="%02d"%(ens)
    # use this file path after I finish ensemble simulations
    # file = pm.ILS_dir()+'/his'+yyyy+"/Qr"+ens_str+"_mean.nc"
    if 'round' in pm.__name__:
        file = pm.ILS_dir()+'/his2024'+"/Qr00"+"_round"+str(pm.val_round(pm.val_mode()))+"_mean.nc"  # use the oneyear simulation to do anomaly
    else:
        file = pm.ILS_dir()+'/his2024'+"/Qr00"+"_mean.nc"  # use the oneyear simulation to do anomaly
    return file
def calc_mean():
    file_path = pm.expout()+'/simmean/'+pm.runname(pm.mode())+'/'
    os.makedirs(file_path,exist_ok=True)
    for num in range(0,pm.ens_mem()+1):
        numch='%03d'%num
        var_mean = read_nc(read_sim_data(num))
        var_flat = var_mean.flatten() 
        var_bin  = var_flat.astype('float32').tobytes()
        # # write .bin
        fileens ="C"+numch
        filebin  = file_path+fileens+'.bin'
        with open(filebin,'wb') as file:
            file.write(var_bin)
