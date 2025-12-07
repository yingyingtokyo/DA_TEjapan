#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#libralies
import os
import numpy as np
import sys
import subprocess
import shutil
#external python codes
import params_case14 as pm

def download_obs(mode=pm.real_mode()):
    obs_code_dir = pm.obs_download_code()
    da_dir = pm.DA_dir()
    if mode == 0:
        return None
    else:
        download_mode, start_yyyy, start_mm, start_dd, start_hh, end_yyyy, end_mm, end_dd, end_hh = pm.download_obs(mode=pm.real_mode())
        validate_mode = pm.val_mode()
        print('obs starting time (JST):',start_yyyy,start_mm,start_dd,start_hh)
        print('obs ending time (JST):',end_yyyy,end_mm,end_dd,end_hh)

        ####### downloading website information ###################
        def download_var(var):
            # subprocess.run(["python3", obs_code_dir + "/1.mlit_html_req.py",start_yyyy,end_yyyy,var]) # 1.only run this code when you are the first time to download data of this year
            subprocess.run(["python3", obs_code_dir + "/2.mlit_wget_wlv_hour.py",start_yyyy, start_mm, start_dd, end_yyyy, end_mm, end_dd, var]) # 2.downloading data
            subprocess.run(["python3", obs_code_dir + "/3.mlit_make_bin_hour.py",start_yyyy, start_mm, start_dd, start_hh, end_yyyy, end_mm, end_dd, end_hh, var]) # 3.convert to bin format
            subprocess.run(["python3", obs_code_dir + "/7.mlit_wlv_hour.py",start_yyyy, start_mm, start_dd, start_hh, end_yyyy, end_mm, end_dd, end_hh]) # 3.convert to hourly format
        
        #download_var('wlv')
        if validate_mode == 2:
            download_var('dis')

        ####### observation allocation ###########################
        ## allocation code is in pm.CaMa_dir()+/map/tej_01min/src_param/allocate_gauge.F90
        ## code for converting observations for generating localization parameters is in ./localization
        localization_dir = da_dir + "/src/localization"
        # subprocess.run(["python3", localization_dir + "/1.mlit_wlv_hour.py",start_yyyy, start_mm, start_dd, start_hh, end_yyyy, end_mm, end_dd, end_hh])  
        ## remake localization parameters
        #if pm.remake_patch() == 1:
        #    os.system(localization_dir + "/3.make_sft.sh "+ str(pm.patch_size())+" " + pm.CaMa_dir()+" " + localization_dir)
        ############################################################

        ###### setting of local patch  ###########################


def prepare_obs():
    dahour_str = '{:02d}'.format(pm.dahour())
    os.system("mkdir -p "+pm.outdir()+"assim_"+pm.runname(pm.mode())+dahour_str)
    if not os.path.islink(pm.outdir()+"assim_"+pm.runname(pm.mode())+dahour_str+"/wlv"):
        os.system("ln -sf "+pm.obs_dir()+" "+pm.outdir()+"assim_"+pm.runname(pm.mode())+dahour_str+"/wlv")
    return 0
####################################
if __name__ == "__main__":
	print ("prepare observations")
	prepare_obs()
