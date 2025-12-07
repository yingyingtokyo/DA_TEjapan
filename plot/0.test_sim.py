#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:06:49 2024
# calculate the range of ensemble simulation

@author: yingyingliu
"""

#%%
import numpy as np
import netCDF4 as nc
import os
import re

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,BoundaryNorm
import matplotlib.cm as cm
import imageio

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.io import shapereader
import os
#%%
ens = 10
savepath = '/work/a06/yingying/Code/plot/exp_tej/'
os.makedirs(savepath, exist_ok=True)

# 100500-102423
day1 = 10
day2 = 18
mon  = 1
var  = 'outflw'
time = (day2-day1)*24
ens_start = 1
mean_num = 10
filepath = '/work/a06/yingying/camada/HydroDA/src/CaMa_out_tej/'

def read_bin(mon,date,enum,chr_type):
    month = "{:02d}".format(mon)
    day   = "{:02d}".format(date)
    ens   = "{:03d}".format(enum)
    fname = filepath+'2019'+month+day+chr_type+ens+'/'+var+'2019.bin'
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(24,1320,1500))
    return data_re
data_test = read_bin(1,1,1,'C')
print(data_test[:,583,1005])
