#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re
import os
from datetime import datetime,timezone,timedelta

########################
#
# parameters list
#
########################
# *************************************************************
def def_params():
    pm_line = 'import params_ils06 as pm\n'
    def replace_as_pm_line(target_file):
        # replace the definition of pm in target files 
        with open(target_file, 'r') as f:
            target_lines = f.readlines()
        new_target_lines = []
        for line in target_lines:
            if 'as pm' in line:
                new_target_lines.append(pm_line)
            else:
                new_target_lines.append(line)
        with open(target_file, 'w') as f:
            f.writelines(new_target_lines)
        print("finish definition of params.py！")

    replace_as_pm_line('main_code.py')
    replace_as_pm_line('prep_patch.py')
    replace_as_pm_line('prep_init.py')
    replace_as_pm_line('prep_runoff.py')
    replace_as_pm_line('prep_obs.py')

# 0. Obs preparation
def mode():
    return 3 
    # parameter to change assimilation mode
    # runoff ensembles will change accordingly.
    # 1: ERA5
    # 2: tej  -- runoff perturbation
    # 3: ils  -- rainfall perturbation
    # 4: topography
    # 5: wid
    # 6: hgt
    # 7: man  -- manning coefficient perturbation
    # 8: msm  -- forecast with rainfall forcing from MSM


def flood_reg():
    # case0 typhoon range: 33.4-41.2; 135.4-142
    ### return [lat_min,lat_max,lon_min,lon_max]
    return [33.4,41.2,135.4,142]
# do 24h forecast
def DA_leading():
    if mode()==8:
        # when it's the first time running this script, it cannot be 'msm'
        return 'on'
    else:
        # when it's the first time running this script, it must be off
        # return "on"
        return "off"

def dahour():
    # time interval for DA
    return 6 

def da_ehour():
    # length of simulation time
    if DA_leading()=='on':
        return 24
    else:
        return dahour() 
        # length of simulation time is the same as the DA interval

def starttime(DAhour=DA_leading()):
    #if DAhour == 'on':
        return (2019,10,7,0) # start date: [year,month,date]
        #return (2019,10,11,0) # start date: [year,month,date]
    #else:
    #    return (2022,1,1)

def endtime(DAhour=DA_leading()):
    #if DAhour == 'on':
        return (2019,10,17,0) # end date: [year,month,date]
        #return (2019,10,12,0) # end date: [year,month,date]
    #else:
    #    return (2022,12,31)  # *note: this date is not included

def real_mode():
    return 0   # turn off real-time mode and don't need to download observations
    # return 1   # turn off real-time mode but download one-year observations
    # return 2   # turn off real-time mode but download one-year observations and remake observation patch
    # return 3 # turn on real-time mode and remake observation patch, for all observations
    # return 4 # turn on real-time mode but not remake observation patch

def val_mode():
    # return 0  # turn off validation mode
    return 1  # turn on validation mode and no need to download discharge for validation
    # return 2  # turn on validation mode but no need to download discharge for validation
    # return 3  # turn on validation mode and download discharge for validation but no need to reselect stations
    # return 4  # turn on validation mode and download discharge for validation
    # return 5  # turn on validation mode and start rolling selection

def val_round(mode=val_mode()):
    if mode==5:
        return 1
    else:
        return 0


def download_obs(mode=real_mode()):
    if (mode == 0) | (mode == 1) | (mode == 2):
        download_mode = "off"
        # start time of observations
        start_yyyy = '{:04d}'.format(starttime(DAhour=DA_leading())[0])
        start_mm = '{:02d}'.format(starttime(DAhour=DA_leading())[1])
        start_dd = '{:02d}'.format(starttime(DAhour=DA_leading())[2])
        start_hh = '{:02d}'.format(starttime(DAhour=DA_leading())[3] + 9)
        # ending time of observations
        end_yyyy = '{:04d}'.format(endtime(DAhour=DA_leading())[0])
        end_mm = '{:02d}'.format(endtime(DAhour=DA_leading())[1])
        end_dd = '{:02d}'.format(endtime(DAhour=DA_leading())[2])
        end_hh = '{:02d}'.format(endtime(DAhour=DA_leading())[3] + 9)
        return download_mode, start_yyyy, start_mm, start_dd, start_hh, end_yyyy, end_mm, end_dd, end_hh   # downloading historical observations
    else:
        download_mode = "on"
        now = datetime.now(timezone.utc)
        start = now - timedelta(days=1)
        start_yyyy = '{:04d}'.format(start.year)
        start_mm = '{:02d}'.format(start.month)
        start_dd  = '{:02d}'.format(start.day)
        start_hh = '{:02d}'.format(start.hour + 9)
        end_yyyy = '{:04d}'.format(now.year)
        end_mm = '{:02d}'.format(now.month)
        end_dd  = '{:02d}'.format(now.day)
        end_hh = '{:02d}'.format(now.hour + 9)        
        return download_mode, start_yyyy, start_mm, start_dd, start_hh, end_yyyy, end_mm, end_dd, end_hh  # downloading real-time observations

def download_dis():
    # return on
    return off

# **************************************************************
# 1. experment type related definitions
def conflag():
    return 2
    # converstion flag for observation converstions 
    #  1 - Directly values 
    #  2 - Anomalies
    #  3 - Normalized values
    #  4 - Log converted values

def mapname():
    return "tej_01min"
    # return "amz_06min"
    # return "glb_15min"
    # realted CaMa-Flood map directory
    # [e.g. : glb_15min, glb_06min, Mkg_06min, etc.]
    # Check 

def remake_patch():
    #return 0  # don't remake patch, just use the default one
    return 1  # remake patch

def patch_size():
    # return 0
    return 100
    # the size of the local patch of LETKF(Local ** EnKF)
    # 0: only 1 pixel (the pixel itself) belongs to its local patch
    # 100: empirical local patch

def DA_dir():
    return "/data42/yingying/HydroDA"
    # directory of HydroDA
    # where src, dat, sat, out exsits

def thersold():
    # return 0.80
    return 0.60
    # return 0.40
    # return 0.20
    # thersold to define the local patch

def patch_dir():
    return "/data42/yingying/HydroDA/obs_patch0"+str(int(thersold()*10))

def patch_filedir():
    os.makedirs("/data42/yingying/HydroDA/obs_patch/obs_patch0"+str(int(thersold()*10)), exist_ok=True)
    return "/data42/yingying/HydroDA/obs_patch/obs_patch0"+str(int(thersold()*10))

def rho():
    # return -1.0
    return 1.00
    # return 1.08
    # -1.0 : adaptive inflation will be used as in Myoshi et al (2011)
    # positive : fixed inflation parameter will be used
    # [E.g. 1.08, 1.10]

def sigma_b():
    # bacground variance of inflation for adaptive inflation Myoshi et al (2011)
    return 0.0400000

def ens_mem(mode=mode()):
    # number of ensemble members
    return 20

def outdir():
    if (val_mode()==0):
        return '/data50/yingying/HydroDA/out/'
    else:
        return '/data50/yingying/HydroDA/out_valid_round'+ str(val_round(val_mode())) +'/'

def expout():
    if (val_mode()==0):
        return '../out'
    else:
        return '../out_valid_round'+ str(val_round(val_mode())) +'/'

# **************************************************************
# 3. Experiment forcing
def timestep():
    return 3600 # outer timestep in seconds

def ILS_dir():
    return '/data42/yingying/ILS/hagibis'

def force_mode(num=mode()):
    if num == 3:
	# rainfall pertubation
        return 1
    if num == 2:
        # runoff pertubation
        return 2
    if num == 7:
        # manning coefficient pertubation
        return 3

# **************************************************************
# 4. Spinup options
def spinup_end_year():
    return 2019

def spinup_end_month():
    return 12

def spinup_end_date():
    return 31

def calibrate():
    return "no"
    # return "yes"
    # return "corrupt"

# **************************************************************
# 5. Runoff forcing 
def run_flag():
    return 0
    # 0: parallel

def runoff_dir():
    return "/data42/yingying/cama/inp"

def runname(num=mode()):
    if num == 1:
        return "ERA5"

    if num == 2:
        return "tej"

    if num == 3:
        return "ils"

    if num == 4:
        return "tpo"

    if num == 5:
        return "wid"

    if num == 6:
        return "hgt"

    if num == 7:
        return "man"

    if num == 8:
        return "ils"

def input(num=mode()):
    if num==1:
        return "ERA5"

    if num==2:
        return "tej"
    # define the runoff data type.

    if num == 3:
        return "ils"

    if num == 4:
        return "tpo"

    if num == 5:
        return "wid"

    if num == 6:
        return "hgt"

    if num == 7:
        return "man"

    if num == 8:
        return "msm"

# **************************************************************
# 6. CaMa-Flood settings
def CaMa_ver():
    # return "CaMa-Flood version 3.9.6"
    return "CaMa-Flood version 4.1.1"

def CaMa_dir():
    return "/data42/yingying/cama/"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders

def CaMa_out(num=mode()):
    if num == 2:
        return "CaMa_out_tej"
    if num == 3:
        return "CaMa_out_ils"
    if num == 4:
        return "CaMa_out_tpo"
    if num == 5:
        return "CaMa_out_wid"
    if num == 6:
        return "CaMa_out_hgt"
    if num == 7:
        return "CaMa_out_man"
    if num == 8:
        return "CaMa_out_ils"
    else:
        return "CaMa_out"

def CaMa_in(num=mode()):
    if num == 2:
        return "CaMa_in_tej"
    if num == 3:
        return "CaMa_in_ils"
    if num == 4:
        return "CaMa_in_tpo"
    if num == 5:
        return "CaMa_in_wid"
    if num == 6:
        return "CaMa_in_hgt"
    if num == 7:
        return "CaMa_in_man"
    if num == 8:
        return "CaMa_in_ils"
    else:
        return "CaMa_in"

def CaMa_out_dir():
    return CaMa_out(num=mode()) +'{:02d}'.format(dahour())

def CaMa_in_dir():
    return CaMa_in(num=mode()) + '{:02d}'.format(dahour()) 

def output_er():
    return 0
    # setting for saving or deleting intermediate files
    # 0 for saving & 1 for deleting
    # those files may be more than 400GB, so erasing is recommended if not necessary

# **************************************************************
# 7. observations settings

def obs_name():
    # return "HydroWeb"
    return "MLIT"

def HydroWeb_dir():
    return "/cluster/data6/menaka/HydroWeb"

def obs_download_code():
    return DA_dir() + "/Empirical_LocalPatch/MLIT"

def obs_dir(mode=real_mode()):
    start_year,start_month,start_date,start_hour=starttime()
    if mode == 0:
        # no need to download observations
        return "/data42/yingying/obs/"+"%04d"%(start_year)+"UST/wlv"
    if (mode == 1) | (mode == 2):
        # download observations from website
        return DA_dir() + "/Empirical_LocalPatch/MLIT/"+"%04d"%(start_year)+'UST/wlv'
    else:
        return DA_dir() + "/Empirical_LocalPatch/MLIT/realtime/wlv"

def alloc_dir():
        return DA_dir() + "/Empirical_LocalPatch/MLIT/alloc"

def plot_dir():
    return DA_dir() + '/plot' 

# **************************************************************
# 8. parallel run settings
def para_nums():
    return 5 
    # setting number of parallels to run CaMa-Flood Model
    # defualt is 6, but may change depending on your system

def cpu_nums():
    return 20 
    # number of cpus used
