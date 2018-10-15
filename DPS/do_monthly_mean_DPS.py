'''
    read DPS model output from the reanalysis run
    and calc the monthly mean and save all the years
    (2002-2016 Jun) in the same file

    ONLY for the Australian region (-55:10N -- 90:180E)

   Author: Eva A Cougnon
   Created: Oct 2017
   Last modification:  

'''

# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import time as time

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_useful as eac

# load data
header_out = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'
outfile = header_out + 'sst_reana_ETKF_mem001_20032017_month_Aus.nc'
outfile_Oz = header_out + 'sst_reana_ETKF_mem001_20032017_daily_Aus.nc'
# define the region of interest
lat_min = -55 
lat_max = 10 
# !!! LONGITUDE in the model!!! from -280 to 80 deg E
lon_min = -270 #90 when deg east from 0 to 360
lon_max = -180 #180 when deg east from 0 to 360
yt_ocean = xr.open_dataset(gname)['yt_ocean']
yt_ocean = yt_ocean.sel(yt_ocean=slice(lat_min,lat_max))
xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
xt_ocean = xt_ocean.sel(xt_ocean=slice(lon_min,lon_max))
#tim = xr.open_dataset(fname)['time']

# VERSION == 1 with the first set of data from Terry -- reanalysis
#   from 2002 - Jun 2016
# VERSION == 2 with the second set of data from Terry -- ETKF

VERSION = 2

if VERSION == 1:
# do not use the year 2002 as it is spinning up
    MinYear = 2003
    MaxYear = 2017
# time vectors
    tim_vec = pd.date_range('2003-01-01','2017-12-30',name='time',freq='M')
    tim_day_vec = pd.date_range('2003-01-01','2017-12-19',name='time',freq='D')
# allocate memory for monthly sst
    X = len(xt_ocean)
    Y = len(yt_ocean)
    nm = len(tim_vec)
    nd = len(tim_day_vec)
    SST_m = xr.Dataset({'sst_m':(('time','yt_ocean','xt_ocean'), \
                                 np.zeros((nm, Y, X)))}, \
                       {'time': tim_vec, 'yt_ocean': yt_ocean, \
                        'xt_ocean': xt_ocean})
    SST_d = xr.Dataset({'sst_d':(('time','yt_ocean','xt_ocean'), \
                                 np.zeros((nd, Y, X)))}, \
                       {'time': tim_day_vec, 'yt_ocean': yt_ocean, \
                        'xt_ocean': xt_ocean})
    k=0
    k_=0
    for i in range(MinYear,MaxYear+1):
        t = time.time()
        print(i)
        fname_i = header + 'sst_reana_' + str(i) + '_daily.nc'

        sst_tmp = xr.open_dataset(fname_i)['sst']
        sst_tmp = sst_tmp.sel(yt_ocean=slice(lat_min,lat_max), \
                              xt_ocean=slice(lon_min,lon_max))

        month_avg_sst = sst_tmp.resample('1MS', dim='time', how='mean')
        if i!=2016:
            SST_m['sst_m'][k:k+12,:,:] = month_avg_sst
        elif i==2016:
#statement for year == 2016! only 6 months!!
            SST_m['sst_m'][k:k+6,:,:] = month_avg_sst
        k = k +12
        tmp = len(sst_tmp[:,0,0])
        SST_d['sst_d'][k_:k_+tmp,:,:] = sst_tmp
        k_ = k_ + tmp
        elapsed = time.time() - t
        print('elapsed time for each year:', elapsed)
elif VERSION == 2:
# daily data for the Australian region
    SST_d = eac.read_netcdfs('/home/ecougnon/data/DPS/reanalysis/ETKF/ocean_daily_SST_enkf9_mem001_20??-20??.nc', \
                             dim='time', lat=yt_ocean, lon=xt_ocean, \
                             transform_func=None)
    SST_d = SST_d.squeeze('st_ocean')
# monthly mean for the Australian region
    SST_m = eac.read_netcdfs('/home/ecougnon/data/DPS/reanalysis/ETKF/ocean_daily_SST_enkf9_mem001_20??-20??.nc', \
                             dim='time', lat=yt_ocean, lon=xt_ocean, \
                             transform_func=lambda ds: ds.resample('1MS', \
                                                               dim='time', \
                                                               how='mean'))
    SST_m = SST_m.squeeze('st_ocean')

    

## saving
SST_m.to_netcdf(outfile)
SST_d.to_netcdf(outfile_Oz)




