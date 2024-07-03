import numpy as np
import os
import numpy 
import math
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import pandas as pd

def treshold_shufling( nfilas,num_surr,dir_surr,id_cluster_1, id_cluster_2):
   
    mxsurr = np.zeros(num_surr)
       
    for columnas in range(num_surr):
        addirsurr = dir_surr+"/surr_"+str(columnas+1)
       
        mx_adir=str(addirsurr+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+'.txt')
              
        dxsr, dysr = np.loadtxt(mx_adir,skiprows=0, usecols=[0,1],unpack=True)   
   
        mxsurr[columnas] = dysr[0]
        print(dysr[0])

    file_txt="significancia.txt"         
    fdata = open(file_txt, 'w')

    try:

        for zz_i in range(len(mxsurr)):
            fdata.write( str(mxsurr[zz_i])+"\n")
                            
    finally:
        fdata.close()
    
    
treshold_shufling( 1,3000,"srExp_sua_2020Mar4w20",393, 427)