'''
    generate and save data to then use in plot_mapX.py
    data correspond to: number of days above the 90th percentile
    or 10th per pixels around Australia as well as the threshold 
    value for each pixel 

	firstly using the whole time series, then using the years 
	  corresponding to a specific climate mode (will be region
	  dependant!!!!)


    Also plot the map of the number of Marine Heat (Cold) Waves
    event using the marine extreme tool box developped by Eric 
    Oliver based on the definition of MHW by Hobbay et al. 2016

   Author: Eva A Cougnon
   Created: Feb 2017
   Last modification: May 2017 (update +ve, -ve limit of ENSO)

'''

# load required modules

import numpy as np
from scipy import stats
import pandas as pd
import xarray as xr

from datetime import date
from calendar import monthrange
import time as time

# Load data
fname = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
tim = xr.open_dataset(fname)['time']
sst_ts = xr.open_dataset(fname)['SSTa']
outfile = '/home/ecougnon/ana/SSTa_map_Aus.nc'

# usefull numbers
MinYear = 1982
MaxYear =2016
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
# time vector for plotting
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,31).toordinal()+1)
# allocate memory
Y = len(lon)
X = len(lat)
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
'''
SST = xr.Dataset({'Tmean':(('lon','lat'),np.zeros((Y, X))), \
                  'T_max':(('lon','lat'),np.zeros((Y, X))), \
                  'T_min':(('lon','lat'),np.zeros((Y, X))), \
                  'Tp90':(('lon','lat'),np.zeros((Y, X))), \
                  'Tp10':(('lon','lat'),np.zeros((Y, X)))}, \
                 {'lat': lat, 'lon': lon})
## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 

file_enso = np.genfromtxt('/home/ecougnon/data/enso/MEI_index.txt', \
                        skip_header=10, skip_footer = 30, delimiter='\t')
str_id = np.nonzero((file_enso[:,0]>(MinYear-1)) \
                    & (file_enso[:,0]<(MinYear+1)))

enso_id = np.zeros(NumYears*12) # the idexes to use for our study (from 1982)
dtime_enso = [None] * (NumYears*12)
enso_daily = np.empty(NumDays)
k=0
l=0
d=0
for yy in np.arange(str_id[0][0],len(file_enso[:,0])):
    for mm in np.arange(1,12+1):
        enso_id[k] = file_enso[yy,mm]
        dtime_enso[k] = date(MinYear+l,mm,1).toordinal()
        enso_daily[d:d+monthrange(MinYear+l,mm)[1]] = enso_id[k]
        k = k + 1
        d = d + monthrange(MinYear+l,mm)[1]
    l = l + 1
enso_p_id = np.nonzero(enso_daily>=0.75)
enso_n_id = np.nonzero(enso_daily<=-0.75)
enso_neu_id = np.nonzero((enso_daily>-0.75) & (enso_daily<0.75))


## basics stats
# SSTs for positive ENSO
'''
sst_enso = np.array(sst_ts) #.copy()
sst_enso[:,:,enso_p_id] = np.nan
sst_enso[:,:,enso_neu_id] = np.nan
'''
# SSTs stats for the entire time series
sst_enso = np.array(sst_ts)

# find p90 and p10 threshold as well as max and min
SST['Tp90'][:,:] = np.apply_along_axis(np.nanpercentile,2,sst_enso[:,:,:],90)
SST['Tp10'][:,:] = np.apply_along_axis(np.nanpercentile,2,sst_enso[:,:,:],10)
SST['T_max'][:,:] = np.apply_along_axis(np.nanmax,2,sst_enso[:,:,:])
SST['T_min'][:,:] = np.apply_along_axis(np.nanmin,2,sst_enso[:,:,:])  
  
## save files
SST.to_netcdf(outfile)



