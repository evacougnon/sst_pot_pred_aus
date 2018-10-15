# load SSTs daily time series and estimate a percentage of the data
# above/below p90/p10 values during El Nino and La Nina yeats
#
# Also estimate a level of significance on the data by ramdomly changing
# the starting point of the ENSO time series (keeping the same statistical
# properties of the ENSO time series) then do the same calc 1000+ times
# and estimate the 95% (or other) significance on the data
#
##############################################################

# load required modules

import numpy as np
from scipy import stats
import pandas as pd
import xarray as xr

from datetime import date
from calendar import monthrange
import time as time

# saving file
sname = '/home/ecougnon/ana/NumDays_ElN_p90p10_Aus_95sign_1000'
sname2 = '/home/ecougnon/ana/NumDays_LaN_p90p10_Aus_95sign_1000_bis'

# Load data
# Daily SSTa
fname = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
#tim = xr.open_dataset(fname)['time']
# number of days for 10% of the total number of days
#tim_10th = len(tim)/10
# SSTa time series
sst_ts = np.array(xr.open_dataset(fname)['SSTa'])
# threshold values (from the daily SSTa time series) at p10 and p90
fname2 = '/home/ecougnon/ana/SSTa_map_Aus.nc'
#SST_X = xr.open_dataset(fname2)
sst_p10 = np.array(xr.open_dataset(fname2)['Tp10'])
sst_p90 = np.array(xr.open_dataset(fname2)['Tp90'])

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
n_p90_LaN -- number of days above p90 corresponding to La Nina year
n_p10_LaN -- number of days below p10 corresponding to La Nina year
'''

Pdays_LaN = xr.Dataset({'n_p90_LaN':(('time','lon','lat'),np.zeros((rtest, Y, X))), \
                        'n_p10_LaN':(('time','lon','lat'),np.zeros((rtest, Y, X)))}, \
                       {'lat': lat, 'lon': lon, 'time':np.arange(0,rtest)})

Pdays_ElN = xr.Dataset({'n_p90_ElN':(('time','lon','lat'),np.zeros((rtest, Y, X))), \
                        'n_p10_ElN':(('time','lon','lat'),np.zeros((rtest, Y, X)))}, \
                       {'lat': lat, 'lon': lon, 'time':np.arange(0,rtest)})

## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 
file_enso = np.genfromtxt('/home/ecougnon/data/enso/MEI_index.txt', \
                        skip_header=10, skip_footer = 30, delimiter='\t')
str_id = np.nonzero((file_enso[:,0]>(MinYear-1)) \
                    & (file_enso[:,0]<(MinYear+1)))

# the indeces to use for our study (from 1982)
enso_id = np.zeros(NumYears*12) 
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
# so it starts at any point of the enso time series and loop back to the 
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

    for ln in range(0,Y):
        for lt in range(0,X):
                    
            Pdays_LaN['n_p90_LaN'][rt,ln,lt] = np.count_nonzero(sst_ma_LaN[ln,lt,:] \
                                                                > sst_p90[ln,lt], \
                                                                axis=2)
            Pdays_LaN['n_p10_LaN'][rt,ln,lt] = np.count_nonzero(sst_ma_LaN[ln,lt,:] \
                                                                < sst_p10[ln,lt], \
                                                                axis=2)
                  
            Pdays_ElN['n_p90_ElN'][rt,ln,lt] = np.count_nonzero(sst_ma_ElN[ln,lt,:] \
                                                                > sst_p90[ln,lt], \
                                                                axis=2)
            Pdays_ElN['n_p10_ElN'][rt,ln,lt] = np.count_nonzero(sst_ma_ElN[ln,lt,:] \
                                                                < sst_p10[ln,lt], \
                                                                axis=2)
            
    del sst_ma_LaN
    del sst_ma_ElN
    elapsed_rt = time.time() - rt_time
    print('elapsed time for each random test:', elapsed_rt)

# low and high (2.5th percentile and 97.5th percentile) value for the 95% 
# level confidence 
rt_low_LaN_p90 = np.apply_along_axis(np.percentile,0,Pdays_LaN['n_p90_LaN'],2.5)
rt_high_LaN_p90 = np.apply_along_axis(np.percentile,0,Pdays_LaN['n_p90_LaN'],97.5)
rt_low_LaN_p10 = np.apply_along_axis(np.percentile,0,Pdays_LaN['n_p10_LaN'],2.5)
rt_high_LaN_p10 = np.apply_along_axis(np.percentile,0,Pdays_LaN['n_p10_LaN'],97.5)

rt_low_ElN_p90 = np.apply_along_axis(np.percentile,0,Pdays_ElN['n_p90_ElN'],2.5)
rt_high_ElN_p90 = np.apply_along_axis(np.percentile,0,Pdays_ElN['n_p90_ElN'],97.5)
rt_low_ElN_p10 = np.apply_along_axis(np.percentile,0,Pdays_ElN['n_p10_ElN'],2.5)
rt_high_ElN_p10 = np.apply_along_axis(np.percentile,0,Pdays_ElN['n_p10_ElN'],97.5)


## save files
np.savez(sname,rt_low_ElN_p90=rt_low_ElN_p90, rt_high_ElN_p90=rt_high_ElN_p90, \
         rt_low_ElN_p10=rt_low_ElN_p10, rt_high_ElN_p10=rt_high_ElN_p10)

np.savez(sname2,rt_low_LaN_p90=rt_low_LaN_p90, rt_high_LaN_p90=rt_high_LaN_p90, \
         rt_low_LaN_p10=rt_low_LaN_p10, rt_high_LaN_p10=rt_high_LaN_p10)



