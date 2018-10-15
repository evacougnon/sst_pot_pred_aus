'''
    checking case study for the reanalyses to check in what 
    extend it reproduces the events (WA 2011 and TAS 2016)
'''

# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import datetime
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

# define indexes for lat lon of one pixel
lat_px = -30 #-42 #-30 # deg N
lon_px = 113 #151 #112 # deg E
lon_px_mdl = lon_px-360 ## offset in model

#figfile = '/home/ecougnon/data/DPS/reanalysis/ts_WA_mdl_obs_detrended.png'

# load data
fname_mdl = '/home/ecougnon/data/DPS/reanalysis/sst_reana_20022016_month_Aus_V2.nc'
lat = xr.open_dataset(fname_mdl)['yt_ocean']
lat = lat.sel(yt_ocean=lat_px, method='nearest')
lon = xr.open_dataset(fname_mdl)['xt_ocean']
lon = lon.sel(xt_ocean=lon_px_mdl, method='nearest')
tim = xr.open_dataset(fname_mdl)['time']
sst_mdl = xr.open_dataset(fname_mdl)['sst_m']
sst_mdl = sst_mdl.sel(yt_ocean=lat,xt_ocean=lon)
'''
valid = ~np.isnan(sst_mdl)
if (valid.any()==True):
    sst_mdl[:] = signal.detrend(sst_mdl[:], axis=0, \
                                type='linear')
elif (valid.all()==False):
    sst_mdl[:] = np.nan
'''

fname_obs = '/home/ecougnon/ana/SST_monthly_extremes_Aus'
# the file is nearly 2gb!!!.... takes a while to load!
data = np.load(fname_obs + '.npz')
lon_obs = data['lon_map']
lat_obs = data['lat_map']
sst_obs_X = data['SST'].item()
# use only the monthly mean
tim_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
tim_str = np.where(tim_vec == np.datetime64(datetime.datetime(2002,1,31)))
tim_end = np.where(tim_vec == np.datetime64(datetime.datetime(2016,6,30)))
sst_obs_tmp = sst_obs_X['TMM'][:,:,tim_str[0][0]:tim_end[0][0]+1]
# swap axes to get time, lat, lon
sst_obs_tmp = np.swapaxes(np.swapaxes(sst_obs_tmp,2,0),1,2)
# change to a xarray to apply the same
X = len(lon_obs)
Y = len(lat_obs)
nm = len(tim)
SST_OBS = xr.Dataset({'sst_obs':(('time','lat','lon'), np.zeros((nm,Y,X)))}, \
                     {'time': tim, 'lat': lat_obs, 'lon': lon_obs})
SST_OBS['sst_obs'][:,:,:]=sst_obs_tmp[:,:,:]
#sst_obs=xr.DataArray(SST_OBS['sst_obs'])
sst_obs = SST_OBS['sst_obs'].sel(lat=lat_px,lon=lon_px, method='nearest')
'''
valid = ~np.isnan(sst_obs)
if (valid.any()==True):
    sst_obs[:] = signal.detrend(sst_obs[:], axis=0, \
                                type='linear')
elif (valid.all()==False):
    sst_obs[:] = np.nan
'''


## plotting
plt.figure(figsize=(13,7))
plt.plot(tim,sst_mdl)
plt.plot(tim,sst_obs)
plt.legend(['mdl','obs'])
plt.title('single point time serries -- near 30S 112E')
plt.grid()
plt.show()
#plt.savefig(figfile, bbox_inches='tight', dpi=300)
 





