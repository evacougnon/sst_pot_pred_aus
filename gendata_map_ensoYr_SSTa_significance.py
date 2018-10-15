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
sst_ts = np.array(xr.open_dataset(fname)['SSTa'])
sname2 = '/home/ecougnon/ana/SSTa_map_nina075Yr_Aus_sig'
sname1 = '/home/ecougnon/ana/SSTa_map_nino075Yr_Aus_sig'

# usefull numbers
MinYear = 1982
MaxYear =2016
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
# number of random test
rtest = 1000
# allocate memory
Y = len(lon)
X = len(lat)
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
'''
SST_LaN = xr.Dataset({'Tp90':(('time','lon','lat'),np.zeros((rtest,Y, X))), \
                      'Tp10':(('time','lon','lat'),np.zeros((rtest,Y, X)))}, \
                     {'lat': lat, 'lon': lon, 'time': np.arange(0,rtest)})
SST_ElN = xr.Dataset({'Tp90':(('time','lon','lat'),np.zeros((rtest,Y, X))), \
                      'Tp10':(('time','lon','lat'),np.zeros((rtest,Y, X)))}, \
                     {'lat': lat, 'lon': lon, 'time': np.arange(0,rtest)})
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
# repeat the ENSO index daily for the corresponding month
        enso_daily[d:d+monthrange(MinYear+l,mm)[1]] = enso_id[k]
        k = k + 1
        d = d + monthrange(MinYear+l,mm)[1]
    l = l + 1
## loop 1000 times to get a random starting index of the ENSO index time series
for rt in range(0,rtest):
    rt_time = time.time()
    print(rt+1, 'of', rtest)
# randomly choose a starting index to 'reshuffle' the enso index time series
# so itr starts at any point of the enso time series and loop back to the 
# start to finally get the same length time series
    ii = np.random.randint(0,len(enso_daily))
    enso_rand = np.append(enso_daily[ii:],enso_daily[:ii])
# find index for positive/negative/neutral phase of enso along the
# time series -- use masked array!!!
    enso_ma_LaN = np.ma.masked_greater(enso_rand,-0.75) # La Nina
    enso_ma_ElN = np.ma.masked_less(enso_rand,0.75) # El Nino
# mask SST time series according to the ENSO index
    sst_ma_LaN=sst_ts[:,:,~enso_ma_LaN.mask].copy()
    sst_ma_ElN=sst_ts[:,:,~enso_ma_ElN.mask].copy()
# find p90 and p10 threshold as well as max and min
    SST_LaN['Tp90'][rt,:,:] = np.apply_along_axis(np.nanpercentile,2, \
                                                  sst_ma_LaN[:,:,:],90)
    SST_LaN['Tp10'][rt,:,:] = np.apply_along_axis(np.nanpercentile,2, \
                                                  sst_ma_LaN[:,:,:],10)
    SST_ElN['Tp90'][rt,:,:] = np.apply_along_axis(np.nanpercentile,2, \
                                                  sst_ma_ElN[:,:,:],90)
    SST_ElN['Tp10'][rt,:,:] = np.apply_along_axis(np.nanpercentile,2, \
                                                  sst_ma_ElN[:,:,:],10)  
    del sst_ma_LaN
    del sst_ma_ElN
    elapsed_rt = time.time() - rt_time
    print('elapsed time for each random test:', elapsed_rt)


# low and high (2.5th percentile and 97.5th percentile) value for the 95% 
# level confidence 
rt_low_LaN_p90 = np.apply_along_axis(np.percentile,0,SST_LaN['Tp90'],2.5)
rt_high_LaN_p90 = np.apply_along_axis(np.percentile,0,SST_LaN['Tp90'],97.5)
rt_low_LaN_p10 = np.apply_along_axis(np.percentile,0,SST_LaN['Tp10'],2.5)
rt_high_LaN_p10 = np.apply_along_axis(np.percentile,0,SST_LaN['Tp10'],97.5)

rt_low_ElN_p90 = np.apply_along_axis(np.percentile,0,SST_ElN['Tp90'],2.5)
rt_high_ElN_p90 = np.apply_along_axis(np.percentile,0,SST_ElN['Tp90'],97.5)
rt_low_ElN_p10 = np.apply_along_axis(np.percentile,0,SST_ElN['Tp10'],2.5)
rt_high_ElN_p10 = np.apply_along_axis(np.percentile,0,SST_ElN['Tp10'],97.5)


## save files
np.savez(sname1,rt_low_ElN_p90=rt_low_ElN_p90, rt_high_ElN_p90=rt_high_ElN_p90, \
         rt_low_ElN_p10=rt_low_ElN_p10, rt_high_ElN_p10=rt_high_ElN_p10)

np.savez(sname2,rt_low_LaN_p90=rt_low_LaN_p90, rt_high_LaN_p90=rt_high_LaN_p90, \
         rt_low_LaN_p10=rt_low_LaN_p10, rt_high_LaN_p10=rt_high_LaN_p10)




