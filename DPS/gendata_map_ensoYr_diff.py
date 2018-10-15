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

from datetime import date
from calendar import monthrange
import time as time

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_ClimateIndex as eac

header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
outfile = header + 'SSTa_stats_nina_map_Aus.nc'
# indicate whether we want +ve, +ve or ignore the NINO3.4 phases
#ENSO = 'pos' # consider only the positive phases
ENSO = 'neg' # consider only the negative phases
#ENSO = 'non' # none of that, consider the entire time series

# define the region
lat_min = -55
lat_max = 10
# !!! LONGITUDE in the model!!! from -280 to 80 deg E
lon_min = 90 - 360
lon_max = 180 - 360
# grid info
gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'
yt_ocean = xr.open_dataset(gname)['yt_ocean']
yt_ocean = yt_ocean.sel(yt_ocean=slice(lat_min,lat_max))
xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
xt_ocean = xt_ocean.sel(xt_ocean=slice(lon_min,lon_max))
# usefull numbers
MinYear = 2003
MaxYear = 2017
NumYear = MaxYear - MinYear+1
# warning finishes in Nov 2017 not Dec
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,11,30).toordinal()+1)

# allocate memory
Y = len(xt_ocean)
X = len(yt_ocean)
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
'''
SSTa_stats = xr.Dataset({'Tmean':(('lon','lat'),np.zeros((X, Y))), \
                         'Tmax':(('lon','lat'),np.zeros((X, Y))), \
                         'Tmin':(('lon','lat'),np.zeros((X, Y))), \
                         'Tp90':(('lon','lat'),np.zeros((X, Y))), \
                         'Tp10':(('lon','lat'),np.zeros((X, Y)))}, \
                        {'lat': yt_ocean, 'lon': xt_ocean})

###################################################
# calc NINO3.4 index from DPS
###################################################
nino34 = eac.nino34_index_dps()
mtime = pd.date_range('2003-01-01','2017-12-01',name='time',freq='M')
# make the index daily to be applied on the daily output
nino34_d= np.nan*np.zeros(len(dtime))
m=0
y=0
d=0
for yy in np.arange(0,NumYear):
    if (yy==NumYear-1):
        for mm in np.arange(1,11+1):
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

############################################################
# load the time series
############################################################
SSTa_d = xr.open_dataset(header + 'ssta_reana_ETKF_mem001_20032017_daily_Aus.nc') \
                        ['dsst_mdl'].sel(time=slice('2003-01-01','2017-11-30'), \
                                         yt_ocean=slice(lat_min,lat_max), \
                                         xt_ocean=slice(lon_min,lon_max))

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




