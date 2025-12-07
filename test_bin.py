import numpy as np
def read_bin2d(fname):
    with open(fname, 'rb') as file:
        data_bin = np.fromfile(file, dtype=np.float32)
        data_re  = np.reshape(data_bin,(2,1320,1500))
    return data_re

assimpath = '../out_valid/CaMa_in_ils01/restart/assim/'
simpath = '../out_valid/CaMa_in_ils01/restart/open/'
filen0 =  'restart2019100720A004.bin'
filen1 =  'restart2019100720C004.bin'
data_re0 = read_bin2d(assimpath+filen0)
data_re1 = read_bin2d(simpath+filen1)
print(data_re0[:,435,1031])
print(data_re1[:,435,1031])

