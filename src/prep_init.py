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

#external python codes
import params_case14 as pm
# ################################################
#
# make folders
# download neal-real-time observations from MLIT website 
# initial inflation parameter rho for assimilation
#
# ################################################
def mkdir(path):
    try:
        os.makedirs(path,exist_ok=True)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
# ################################################
def initial(): #used
    # program for initialization
    # creating output folders
    dahour = pm.dahour()
    dahour_str = '{:02d}'.format(dahour)
    mkdir(pm.outdir()+"/logout")
    if not os.path.exists(pm.outdir()+"/CaMa_out_"+pm.runname(pm.mode())):
        os.symlink("/data50/yingying/HydroDA/out_valid/CaMa_out_"+pm.runname(pm.mode()),pm.outdir()+"/CaMa_out_"+pm.runname(pm.mode()))
    mkdir(pm.outdir()+"/CaMa_out_"+pm.runname(pm.mode()))
    mkdir(pm.outdir()+"/CaMa_out_"+pm.runname(pm.mode())+dahour_str)
    mkdir(pm.outdir()+"/CaMa_in_"+pm.runname(pm.mode())+dahour_str)
    mkdir(pm.outdir()+"/CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart")
    mkdir(pm.outdir()+"/CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart/assim")
    mkdir(pm.outdir()+"/CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart/open")
    mkdir(pm.outdir()+"/CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart/true")
    src_fold = pm.runoff_dir() 
    link_fold = pm.outdir()+'/CaMa_in_'+pm.runname(pm.mode())+dahour_str+'/inp'
    if os.path.islink(link_fold) or os.path.exists(link_fold):
            os.remove(link_fold)  # Remove the existing link
    os.symlink(src_fold,link_fold)

    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str)
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m")
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m/assim")
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m/open")
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m/true")

    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa")    
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/assim")  
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/open")
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/meanA")   
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/meanC")  

    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/nonassim")
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/nonassim/open")
    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/nonassim/assim")

    mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/rest_true")
    
    for var in ["rivdph","sfcelv","outflw","storge"]:
        mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/"+var)
        mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/"+var+"/assim")
        mkdir(pm.outdir()+"/assim_"+pm.runname(pm.mode())+dahour_str+"/"+var+"/open")

    # error output folder
    mkdir("err_out")

    # inflation parameter
    mkdir("inflation")

    os.system("touch "+ pm.outdir()  +"/assim_"+pm.runname(pm.mode())+dahour_str+"/__init__.py")
    mkdir("logout")
###########################
def make_initial_infl():
    dahour = pm.dahour()
    dahour_str = '{:02d}'.format(dahour)
    parm_infl=np.ones((1320,1500),dtype=np.float32)
    start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
    yyyy='%04d' % (start_year)
    mm='%02d' % (start_month)
    dd='%02d' % (start_date)
    hh='%02d' % (start_hour)
    os.makedirs(pm.outdir()+'/CaMa_out_'+pm.runname(pm.mode())+dahour_str+"/inflation/",exist_ok=True)
    parm_infl.tofile(pm.outdir()+'/CaMa_out_'+pm.runname(pm.mode())+dahour_str+"/inflation/parm_infl"+yyyy+mm+dd+hh+".bin")
    return 0

