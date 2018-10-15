'''
    checking case study for the reanalyses to check in what 
    extend it reproduces the events (WA 2011 and TAS 2016)

    using the daily output

'''

# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import datetime
from scipy import io # load matlab file
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

# define indexes for lat lon of one pixel
lat_px_min = -44 #-45 #-32 #-45 # deg N
lat_px_max = -42 # -37 #-28 #-37
lon_px_min = 44 #147 #112 #147 # deg E
lon_px_max = 146 #155 #115 #155 
lon_px_min_mdl = lon_px_min-360 ## offset in model
lon_px_max_mdl = lon_px_max-360 ## offset in model

# time vector 
tim_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='D')
tim_str = np.where(tim_vec == np.datetime64(datetime.datetime(1982,1,1)))
tim_end = np.where(tim_vec == np.datetime64(datetime.datetime(2016,12,31)))

#figfile = '/home/ecougnon/data/DPS/reanalysis/ETKF/ts_TAS_area_mdl_obs_SSTadaily.png'

# load data
fname_mdl = '/home/ecougnon/data/DPS/reanalysis/ETKF/ssta_reana_ETKF_mem001_20032017_daily_Aus.nc'
lat_mdl = xr.open_dataset(fname_mdl)['yt_ocean']
lat_mdl = lat_mdl.sel(yt_ocean=slice(lat_px_min,lat_px_max)) #, method='nearest')
lon_mdl = xr.open_dataset(fname_mdl)['xt_ocean']
lon_mdl = lon_mdl.sel(xt_ocean=slice(lon_px_min_mdl,lon_px_max_mdl)) #, method='nearest')
tim_mdl = xr.open_dataset(fname_mdl)['time']
tim_mdl = tim_mdl.sel(time=slice('2003-01-01', '2016-12-31'))
sst_mdl = xr.open_dataset(fname_mdl)['dsst_mdl'] #['temp'] #['dsst_mdl'] #['sst_d']
sst_mdl = sst_mdl.sel(yt_ocean=lat_mdl,xt_ocean=lon_mdl, time=tim_mdl)
sst_mdl = np.nanmean(np.nanmean(sst_mdl,1),1)
'''
valid = ~np.isnan(sst_mdl)
if (valid.any()==True):
    sst_mdl[:] = signal.detrend(sst_mdl[:], axis=0, \
                                type='linear')
elif (valid.all()==False):
    sst_mdl[:] = np.nan
'''
## read obs
## whether or not using the SSTa_daily_Aus.nc
# check = 0 -- pick a file from the original SST.mat files
# check = 1 -- read a nc file 
check = 0

if (check == 0):
    header_obs = '/home/data/sst/noaa_oi_v2/avhrr/'
    matobj_tmp = io.loadmat(header_obs + \
                            'timeseries/avhrr-only-v2.ts.0001.mat')
    lat_tmp = matobj_tmp['lat']
    lon_tmp = matobj_tmp['lon']
    res = 0.25 # resolution of the data -- quarter degree
# find the closest index from the lat/lon defined above
    lat_px_min_id = np.nonzero((lat_tmp>(lat_px_min-res)) & \
                               (lat_tmp<(lat_px_min)))
    lat_px_max_id = np.nonzero((lat_tmp>(lat_px_max-res)) & \
                               (lat_tmp<(lat_px_max)))
    lon_px_min_id = np.nonzero((lon_tmp>(lon_px_min-res)) & \
                               (lon_tmp<(lon_px_min)))
    lon_px_max_id = np.nonzero((lon_tmp>(lon_px_max-res)) & \
                               (lon_tmp<(lon_px_max)))
# load time series from the given point
    sst_obs = np.empty((len(np.arange(lon_px_min_id[0][0], \
                        lon_px_max_id[0][0]+1)), \
                        len(np.arange(lat_px_min_id[0][0], \
                        lat_px_max_id[0][0]+1)), len(tim_vec)))
    sst_obs.fill(np.nan)
    k=0
    for i in range(lon_px_min_id[0][0],lon_px_max_id[0][0]+1):
        print(i)
        matobj = io.loadmat(header_obs + 'timeseries/avhrr-only-v2.ts.' \
                            + str(i).zfill(4) + '.mat')
        sst_obs[k,:,:] = matobj['sst_ts'][lat_px_min_id[0][0]:lat_px_max_id[0][0]+1, \
                                          tim_str[0][0]:tim_end[0][0]+1]
        k = k+1
    sst_obs = np.nanmean(np.nanmean(sst_obs,0),0)
    '''
    valid = ~np.isnan(sst_obs)
    if (valid.any()==True):
        sst_obs[:] = signal.detrend(sst_obs[:], axis=0, \
                                    type='linear')
    elif (valid.all()==False):
        sst_obs[:] = np.nan
    '''
elif (check == 1):
    import pandas as pd
    import xarray as xr
    fname = '/home/ecougnon/ana/SSTa_daily_Aus_20032016Dec.nc'
    lat_obs = xr.open_dataset(fname)['lat']
#    lat_obs = lat_obs.sel(lat=slice(lat_px_min,lat_px_max))
    lat_obs = lat_obs.sel(lat=lat_mdl, method='nearest')
    lon_obs = xr.open_dataset(fname)['lon']
#    lon_obs = lon_obs.sel(lon=slice(lon_px_min,lon_px_max))
    lon_obs = lon_obs.sel(lon=lon_mdl+360, method='nearest')
    tim_obs = xr.open_dataset(fname)['time']
    tim_obs = tim_obs.sel(time=slice('2003-01-01','2016-12-31'))
    sst_obs = xr.open_dataset(fname)['SSTa']
    sst_obs = sst_obs.sel(time=tim_obs, lat=lat_obs,lon=lon_obs)
    sst_obs = np.nanmean(np.nanmean(sst_obs,0),0)

# pearson correlation
PearC, tmp = st.pearsonr(sst_mdl, sst_obs) 

## plotting
plt.figure(figsize=(13,7))
ax = plt.subplot(111)
#plt.plot(tim_mdl,sst_mdl)
plt.plot(tim_vec,sst_obs)
#plt.legend(['mdl','obs'])
#plt.title('SSTa TAS area time series -- 37-45S 147-155E')
plt.title('SSTa TAS area time series -- 42-44S 144-146E')
plt.grid()
#plt.text(0.3, 0.1, 'Pearson Correlation coefficient:' + str(round(PearC,3)), \
#         ha='center', va='center', transform=ax.transAxes, \
#         fontsize=14)

#plt.savefig(figfile, bbox_inches='tight', dpi=300)

plt.show() 





