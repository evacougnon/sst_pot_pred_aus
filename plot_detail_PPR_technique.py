'''
    plot time series to show how the PPR method works

    Author: Eva C.
    Created: Jul 2020
    Last Modif: 

'''
#########################
# load required modules
#########################

import numpy as np
import xarray as xr
from scipy import signal
from scipy import io
from scipy import stats
import pandas as pd
import time as time
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from matplotlib.dates import MonthLocator, YearLocator

fname = '/home/ecougnon/Documents/PotPred_paper/data/SSTa_daily_Aus.nc'
figfile = '/home/ecougnon/Documents/PotPred_paper/figures/PPR_example.png'

#########################################
# define indexes for lat lon of one pixel
#########################################
lat_px = -30 # deg N
lon_px = 110 #150 # deg E

tau = 365 # chunk size for PPR calc.

tim = pd.date_range('1982-01-01','2016-12-31')
tim_yearly = pd.date_range('1982-07-01','2016-07-01', freq='AS-JUL')
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
NumDays = len(tim)

ts = xr.open_dataset(fname)['SSTa'].sel(time=tim,lat=lat_px,lon=lon_px, method="nearest")

# define the start/end indexes for each chunk
str_id = range(0,NumDays,tau)
end_id = range(tau-1,NumDays+1,tau)
NumChunk = len(str_id)

###########################################################
# calc PotPred_vStorchZwiers for the chosen timeseries
##########################################################
# checking that the time series is not empty
if ((np.isnan(np.min(ts)) == False) & \
    (np.isnan(np.max(ts)) == False)):
# chunk mean allocation
    t_j=np.empty(NumChunk)
    t_j.fill(np.nan)
# periodogram of each chunk
    pxx_j = np.empty([NumChunk,tau])#np.empty([NumChunk,int(tau/2)+1])
    pxx_j.fill(np.nan)# power spectrum of each chunk
    f_j = np.empty([NumChunk,tau])#np.empty([NumChunk,int(tau/2)+1])
    f_j.fill(np.nan)# corresponding frequencies
    for c in range(0,NumChunk): # run through each chunk
# chunk meam
        t_j[c] = np.nanmean(ts[int(str_id[c]):int(end_id[c])+1])

        f_j[c,:],pxx_j[c,:]=signal.periodogram(ts[int(str_id[c]) \
                                                  :int(end_id[c])+1],detrend=False, \
                                               return_onesided=False, scaling='density')
# mean of all the chunks
    t_jall = np.nanmean(t_j)
# total inter chunk variance
    tmp=np.zeros(NumChunk)
    for c in range(0,NumChunk):
        tmp[c] = (t_j[c]-t_jall)**2
    Var_interC = 1/(NumChunk-1) * np.nansum(tmp)
# inter chunk variance of the fast varying process
# average periodograms in time -- estimate the spectrum for all the chunks
    pxx_jall=np.mean(pxx_j,axis=0)

    idx = np.argmin(np.abs(f_j[0][:] - (1/tau)))
# variance of the fast varying component
    if idx!=0:
        Var_noise = (1/tau) * pxx_jall[idx]*2   # when double sided
    elif idx==0:
        Var_noise = (1/tau) * pxx_jall[idx]
# variance ratio
    p_vSZ = Var_interC / Var_noise
# testing the null hypothesis to find the significance
    F90 = stats.f.ppf(0.90,tau-1,2*NumChunk)
    F95 = stats.f.ppf(0.95,tau-1,2*NumChunk)
# PPR to plot (normalised between [0:1]: (Var_interC - Var_noise)/Var_interC
p = (Var_interC - Var_noise)/Var_interC 
#########################################################################
# PLOTTING (ts, p, Var_interC (total inter chunk variance), 
#           Var_noise (variance of the fast varying component), 
#           chunk limits, chunk meam (mean of all the chunks)
#########################################################################
fig, ax = plt.subplots(figsize=(14, 6))
ts.plot.line(color = 'k', alpha = 0.4, label = 'SSTa')
plt.hlines(p,'1982-01','2017-01', color = 'r', linestyle='solid', label = 'P')
plt.plot(tim_yearly, t_j, marker='o', color='k', label = 'Mean of each chunks')
plt.hlines(Var_interC,'1982-01','2017-01', color = 'b', linestyle='--', label = 'Variance of the chunk-mean')
plt.hlines(Var_noise,'1982-01','2017-01', color = 'c', linestyle='--', label = 'unpredictable noise component')
ax.set_xlim('1982-01','2017-01')
ax.set_ylim(-4, 4)
yloc = YearLocator()
ax.xaxis.set_major_locator(yloc)
#plt.xticks(pd.date_range('1982-01-01','2016-01-01', freq='2AS'))
ax.grid(True)
ax.legend()

plt.savefig(figfile,bbox_inches='tight', format='png', dpi=300)
