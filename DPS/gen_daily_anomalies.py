'''

	Generate SSTa anomalies from the DPS reanalysis simulation

    Created: Nov 2017
    Author: Eva C.
    Last modification: Jun 2018 -- adapted to calculated the SSTa
			on the forecast using the climatology 
			from enkf-9 (reanalysis)

'''
# Load required modules

import numpy as np
from scipy import io # read/load mat files
from scipy import signal  # detrend
import xarray as xr

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo


outfile = '/home/ecougnon/data/DPS/forecast/ssta_ens_climenkf9_daily_Aus.nc'

# load data
fname_mdl = '/home/ecougnon/data/DPS/reanalysis/ETKF/sst_reana_ETKF_mem001_20032017_daily_Aus.nc'
yt_ocean = xr.open_dataset(fname_mdl)['yt_ocean']
xt_ocean = xr.open_dataset(fname_mdl)['xt_ocean']
tim = xr.open_dataset(fname_mdl)['time']
tim = tim.sel(time=slice('2003-01-01','2017-12-19'))
sst_mdl = xr.open_dataset(fname_mdl)['temp']
sst_mdl = sst_mdl.sel(time=tim)
fname_for = '/home/ecougnon/data/DPS/forecast/sst_daily_Aus_ensemble.nc'
sst_for = xr.open_dataset(fname_for)['sst']. \
             sel(yt_ocean=yt_ocean,xt_ocean=xt_ocean, \
                 time=slice('2003-01-01','2010-10-25'))
time_for =  sst_for.time

#
# allocate memory for the main variable
#
X = len(xt_ocean)
Y = len(yt_ocean)
nd = len(tim)
dsst_mdl = xr.Dataset({'dsst_mdl':(('time','yt_ocean','xt_ocean'), \
                                   np.zeros((nd, Y, X)))}, \
                      {'time': tim, 'yt_ocean': yt_ocean, 'xt_ocean': xt_ocean})
ssta_for = xr.Dataset({'ssta_for':(('time','yt_ocean','xt_ocean'), \
                                   np.zeros((len(time_for), Y, X)))}, \
                      {'time': time_for, 'yt_ocean': yt_ocean, \
                       'xt_ocean': xt_ocean})
#############################################################
# deseason
# remove the seasonal signal estimated from the daily output
##############################################################
dsea, sea, beta = np.apply_along_axis(eo.deseason_harmonic,0, \
                                      sst_mdl[:,:,:],4,365)
'''
# detrend the time series
#dsst_mdl = np.empty((len(tim), len(lat), len(lon)))
#dsst_mdl.fill(np.nan)
for ii in range(0,Y):
    for jj in range(0,X):
        valid = ~np.isnan(sst_mdl[:,ii,jj])
        if (valid.any()==True):
            dsst_mdl['dsst_mdl'][:,ii,jj] = signal.detrend(np.squeeze(dsea[ii,jj]))
        elif (valid.all()==False):
            dsst_mdl['dsst_mdl'][:,ii,jj] = np.nan
'''
for ii in range(0,Y):
    for jj in range(0,X):
        for tt in range(0,len(time_for)):
            ssta_for['ssta_for'][tt,ii,jj] = sst_for[tt,ii,jj] - sea[ii,jj][tt,0]

#####################
## saving
#####################
#dsst_mdl.to_netcdf(outfile)
ssta_for.to_netcdf(outfile)












