'''

	filter spatially the daily SST satellite based observations
	from 1/4 degree AVHRR OISST to the specified window span 
	(e.g. 1 degree box car window)

	Author: Eva C.
	Created: May 2018
	Last modified:

'''

###################################
# load required modules
###################################

import numpy as np
import xarray as xr
import pandas as pd
import time as time

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

###################################
# define indexes for lat lon of one pixel
###################################
lat_px_min = -5 -1 #-55 - 1  # deg N
lat_px_max = 5 + 1 #10 + 1 
lon_px_min = 160 -1 #90 - 1 # deg E
lon_px_max = 210 + 1 #180 + 1 
# 150:285E to cover the entire Equatorial Pacific Ocean for any NINO index
# however running out of mem!!
header = '/home/data/sst/noaa_oi_v2/new_avhrr/'
outfile = '/home/ecougnon/ana/SST_smooth_1deg_NINO4.nc'

t_read = time.time()

'''
ds = xr.open_mfdataset(header + 'sst.day.mean.????.v2.nc')['sst']. \
        sel(time=slice('1982-01-01','2017-12-31'), \
            lat=slice(lat_px_min,lat_px_max), \
            lon=slice(lon_px_min,lon_px_max))
ds = xr.concat(ds, 'time').squeeze()
ds.to_netcdf('/home/ecougnon/ana/SST_NINO4region_tmp.nc')
elapsed_1 = time.time() - t_read
print('elapsed time: ', \
      elapsed_1)

'''

ds_smooth1 = xr.open_dataset('/home/ecougnon/ana/SST_NINO4region_tmp.nc')['sst']. \
                sel(time=slice('1982-01-01','2001-12-31')). \
                rolling(center=True, lat=5).mean(). \
                rolling(center=True, lon=5).mean()
ds_smooth2 = xr.open_dataset('/home/ecougnon/ana/SST_NINO4region_tmp.nc')['sst']. \
                sel(time=slice('2002-01-01','2017-12-31')). \
                rolling(center=True, lat=5).mean(). \
                rolling(center=True, lon=5).mean()
ds_smooth = xr.concat([ds_smooth1, ds_smooth2], dim ='time')
ds_smooth.to_netcdf(outfile)


elapsed_1 = time.time() - t_read
print('elapsed time: ', \
      elapsed_1)

#'''
'''
testing phase on one file only!
header = '/home/data/sst/noaa_oi_v2/new_avhrr/'
ds_oisst_d = xr.open_dataset(header + 'sst.day.mean.2016.v2.nc')['sst']. \
                sel(lat=slice(lat_px_min,lat_px_max), \
                    lon=slice(lon_px_min,lon_px_max)). \
                rolling(center=True, lat=5).mean(). \
                rolling(center=True, lon=5).mean()
'''


