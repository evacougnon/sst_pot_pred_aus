'''
    generate and save data to then use in plot_mapX.py
    data correspond to: number of days above the 90th percentile
    or 10th per pixels around Australia as well as the threshold 
    value for each pixel 

        firstly using the whole time series, then using the years 
          corresponding to a specific climate mode (will be region
          dependant!!!!)


	based on gendata_map_ensoYr_SSTa.py developed for the obs

   Author: Eva A Cougnon
   Created: Apr 2018
   Last modification: 

'''

# load required modules

import numpy as np
from scipy import stats
import pandas as pd
import xarray as xr
from scipy import signal  # detrend

from datetime import date
from calendar import monthrange
import time as time

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_ClimateIndex as eac

# indicate whether we want +ve, +ve or ignore the NINO3.4 phases
#ENSO = 'pos' # consider only the positive phases
ENSO = 'pos' # consider only the negative phases

header = '/home/ecougnon/ana/'
outfile = header + 'SSTa_stats_nina_map_Aus_filter.nc'
if ENSO == 'pos':
    outfile = header + 'SSTa_stats_nino_map_Aus_filter.nc'
elif ENSO == 'neg': 
    outfile = header + 'SSTa_stats_nina_map_Aus_filter.nc'

#ENSO = 'non' # none of that, consider the entire time series

# define the region
lat_min = -55
lat_max = 10
lon_min = 90 
lon_max = 180 
# usefull numbers
MinYear = 2003
MaxYear = 2017
NumYear = MaxYear - MinYear+1
# warning finishes in Nov 2017 not Dec
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,10,31).toordinal()+1)
NumDays = len(dtime)
# Daily SSTa
SSTa_d = xr.open_dataset(header + 'SSTa_daily_filter_Aus_20032017.nc')['SSTa']
lat_map = SSTa_d.lat
lon_map = SSTa_d.lon

# allocate memory
Y = len(lat_map)
X = len(lon_map)
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
'''
SSTa_stats = xr.Dataset({'Tmean':(('lat','lon'),np.zeros((Y, X))), \
                         'Tmax':(('lat','lon'),np.zeros((Y, X))), \
                         'Tmin':(('lat','lon'),np.zeros((Y, X))), \
                         'Tp90':(('lat','lon'),np.zeros((Y, X))), \
                         'Tp10':(('lat','lon'),np.zeros((Y, X)))}, \
                        {'lat': lat_map, 'lon': lon_map})

########################################
## CLIMATE MODE
########################################
######################################
# NINO34 index from the smoothed (1/4 degree to 1degree)
# NOAA OISSTV2 daily dataset
#########################################
ds_oisst = xr.open_dataset('/home/ecougnon/ana/SST_smooth_1deg_NINO34region.nc')['__xarray_dataarray_variable__']. \
              sel(lat=slice(-5,5), \
                  lon=slice(190,240), \
                  time=slice('2003-01-01','2017-10-31')). \
              resample('1MS', dim='time', how='mean')
clim_oisst = ds_oisst.groupby('time.month').mean('time')
ssta_oisst = ds_oisst.groupby('time.month') - clim_oisst
dssta_oisst = np.apply_along_axis(signal.detrend,0,ssta_oisst, type='linear')
nino34 = np.mean(dssta_oisst,axis=(1,2))
mtime = pd.date_range('2003-01-01','2017-10-31',name='time',freq='M')

# warning finishes in Nov 2017 not Dec
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,10,31).toordinal()+1)
# make the index daily to be applied on the daily output
nino34_d= np.nan*np.zeros(len(dtime))
m=0
y=0
d=0
for yy in np.arange(0,NumYear):
    if (yy==NumYear-1):
        for mm in np.arange(1,10+1):
            nino34_d[d:d+monthrange(MinYear+y,mm)[1]] = nino34[m]
            m = m + 1
            d = d + monthrange(MinYear+y,mm)[1]
    else:
        for mm in np.arange(1,12+1):
            nino34_d[d:d+monthrange(MinYear+y,mm)[1]] = nino34[m]
            m = m + 1
            d = d + monthrange(MinYear+y,mm)[1]
    y = y + 1
# save indexes where enso is +ve (p) and -ve (n)
nino34_p_id = np.nonzero(nino34_d>=0.4)
nino34_n_id = np.nonzero(nino34_d<=-0.4)
nino34_neu_id = np.nonzero((nino34_d>-0.4) & (nino34_d<0.4))


######################################################
# basic stats
######################################################
if ENSO == 'pos':
    sst_enso = np.array(SSTa_d) 
    sst_enso[nino34_n_id,:,:] = np.nan
    sst_enso[nino34_neu_id,:,:] = np.nan
elif ENSO == 'neg':
    sst_enso = np.array(SSTa_d)
    sst_enso[nino34_p_id,:,:] = np.nan
    sst_enso[nino34_neu_id,:,:] = np.nan
elif ENSO == 'non':
# SSTs stats for the entire time series
    sst_enso = np.array(SSTa_d)

# find p90 and p10 threshold as well as max and min
SSTa_stats['Tp90'][:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                              sst_enso[:,:,:],90)
SSTa_stats['Tp10'][:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                              sst_enso[:,:,:],10)
SSTa_stats['Tmax'][:,:] = np.apply_along_axis(np.nanmax,0,sst_enso[:,:,:])
SSTa_stats['Tmin'][:,:] = np.apply_along_axis(np.nanmin,0,sst_enso[:,:,:])
SSTa_stats['Tmean'][:,:] = np.apply_along_axis(np.nanmean,0,sst_enso[:,:,:])

## save files
SSTa_stats.to_netcdf(outfile)




