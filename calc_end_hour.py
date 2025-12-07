#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import sys
from sys import *
from datetime import datetime, timedelta

argvs = sys.argv
yyyy=int(argvs[1])
mm=int(argvs[2])
dd=int(argvs[3])
hh=int(argvs[4])
ehour=int(argvs[5])
obj=argvs[6]

yyyymmddhh=f"{yyyy:04d}{mm:02d}{dd:02d}{hh:02d}"
ctime = datetime.strptime(yyyymmddhh, "%Y%m%d%H")
ntime = ctime + timedelta(hours=ehour)
time_new = ntime.strftime("%Y%m%d%H")

if obj=="year":
    print(ntime.strftime('%Y'))
if obj=="month":
    print(ntime.strftime('%m'))
if  obj=="date":
    print(ntime.strftime('%d'))
if  obj=="hour":
    print(ntime.strftime('%H'))
