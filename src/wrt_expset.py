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
from multiprocessing import sharedctypes
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
import params as pm
###########################
def write_text():
    fdir = pm.outdir() + pm.CaMa_out_dir()
    os.makedirs(fdir, exist_ok=True)
    with open(fdir + "/experimental_settings.log", "w") as f:
        f.write("# Experimental Settings\n")
        f.write("#======================================================\n")
        # Experimental Settings
        f.write("# Experiment Mode: "+"%d"%(pm.mode())+"\n")
        # Time domain for analysis
        f.write("# Start Date: %04d-%02d-%02d-%02d\n"%(pm.starttime()))
        f.write("# End Date: %04d-%02d-%02d-%02d\n"%(pm.endtime()))
        # f.write("# Parameter Corruption: "+corruption(pm.corrupt())+"\n")
        # Assimilation Settings
        f.write("# Created at : "+str(datetime.datetime.now()))
    return 0
###########################
if __name__=="__main__":
    write_text()
