'''

    load from the computed single spectral analysis on PCs time series
    determined after applying and EOF analysis (using cdo) 
    on the SSTa on filtered and unfiltered data, then compare
    the RCs

    Following the idea in Monselesan et al 2015 and O'Kane et al 2016

    With the help of -- copied from!:
    https://dept.atmos.ucla.edu/tcd/ssa-tutorial-matlab

    Author: Eva C.
    Created: Oct 2018
    Last Modif:

'''

# libriries
import numpy as np
import xarray as xr
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays
from scipy import signal
from scipy import stats
import pandas as pd
import csv
import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt


########################
# load data
########################
# EOF on the yearly mean (filtered) and 1x1 deg resolution grid
fname_pcf = '../../ana/PotPred/vSZ/SSA/All_SSA_RCs_YearlyOnMonthly.nc'
tim_vec_f = xr.open_dataset(fname_pcf)['time']


# EOF on the unfiltered (daily) SSTa and (1/4)x(1/4) deg resolution grid
fname_pcu = '../../ana/PotPred/vSZ/SSA/All_SSA_RCs_Daily.nc'
tim = pd.date_range('1982-01-01','2016-12-31',name='time',freq='D')
# remove the last day of the year when leap year
# Need a function!!!....
tim = tim[tim !='1984-12-31']
tim = tim[tim !='1988-12-31']
tim = tim[tim !='1992-12-31']
tim = tim[tim !='1996-12-31']
tim = tim[tim !='2000-12-31']
tim = tim[tim !='2004-12-31']
tim = tim[tim !='2008-12-31']
tim = tim[tim !='2012-12-31']
tim = tim[tim !='2016-12-31']

#################################
# which mode to compare with
################################
# correct the time vector to appky the monthly mean
RC_allPCs = xr.Dataset({'RC_allPCs':(('slow_RC','time','modes'), \
                                     xr.open_dataset(fname_pcu) \
                                     ['RC_allPCs'][0:3,:,0:6])}, \
                       coords = {'time':tim, 'slow_RC':np.arange(0,3), \
                                 'modes':np.arange(0,6)})

var_PC_u = RC_allPCs['RC_allPCs'].resample(time='1MS').mean('time')
var_PC_f = xr.open_dataset(fname_pcf)['RC_allPCs'][0:3,:,0:6] 

diff_RCs = var_PC_u - var_PC_f

std_diff_RCs = np.std(diff_RCs,axis=1)


#################################
# plotting
##################################
ax=plt.figure(figsize=(11,3))
plt.plot(var_PC_u[0,:,0]) 
plt.plot(var_PC_f[0,:,0])
plt.legend(['unfiltered','filtered'])
plt.grid()


ax=plt.figure(figsize=(11,3))
k=0
for mm in range(0,6):
    plt.plot(np.arange(k,k+3),std_diff_RCs[:,mm],'*')
    k = k + 3
plt.grid()



plt.show(block=False)



