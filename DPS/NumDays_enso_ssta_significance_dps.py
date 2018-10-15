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

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_ClimateIndex as eac

header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
# saving file
sname = header + 'NumDays_ElN_p90p10_Aus_95sign_2500_dps'
sname2 = header + 'NumDays_LaN_p90p10_Aus_95sign_2500_dps'

################################
# define the region
################################
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
NumDays = len(dtime)
# number of days for 10% of the total number of days
tim_10th = len(dtime)/10
# Load data
# Daily SSTa
SSTa_d = xr.open_dataset(header + 'ssta_reana_ETKF_mem001_20032017_daily_Aus.nc') \
                        ['dsst_mdl'].sel(time=slice('2003-01-01','2017-11-30'), \
                                         yt_ocean=slice(lat_min,lat_max), \
                                         xt_ocean=slice(lon_min,lon_max))
# threshold values (from the daily SSTa time series)
fname2 = header + 'SSTa_stats_map_Aus.nc'
#SST_X = xr.open_dataset(fname2)
sst_p10 = np.array(xr.open_dataset(fname2)['Tp10'])
sst_p90 = np.array(xr.open_dataset(fname2)['Tp90'])
# number of random test
rtest = 2500

####################################
# allocate memory
####################################
Y = len(xt_ocean)
X = len(yt_ocean)
'''
n_p90_LaN -- number of days above p90 corresponding to La Nina year
n_p10_LaN -- number of days below p10 corresponding to La Nina year
'''

Pdays_LaN = xr.Dataset({'n_p90_LaN':(('time','lon','lat'),np.zeros((rtest,X,Y))), \
                        'n_p10_LaN':(('time','lon','lat'),np.zeros((rtest,X,Y)))}, \
                       {'lat': yt_ocean, 'lon': xt_ocean, 'time':np.arange(0,rtest)})

Pdays_ElN = xr.Dataset({'n_p90_ElN':(('time','lon','lat'),np.zeros((rtest,X,Y))), \
                        'n_p10_ElN':(('time','lon','lat'),np.zeros((rtest,X,Y)))}, \
                       {'lat': yt_ocean, 'lon': xt_ocean, 'time':np.arange(0,rtest)})

########################################
## CLIMATE MODE
########################################
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

##############################################
# randomising the time series rtest times
#############################################
for rt in range(0,rtest):
    rt_time = time.time()
    print(rt+1, 'of', rtest)
# randomly choose a starting index to 'reshuffle' the enso index time series
# so it starts at any point of the enso time series and loop back to the 
# start to finally get the same length time series
    ii = np.random.randint(0,len(nino34_d))
    enso_rand = np.append(nino34_d[ii:],nino34_d[:ii]) 
# find index for positive/negative/neutral phase of enso along the
# time series -- use masked array!!!
    enso_ma_LaN = np.ma.masked_greater(enso_rand,-0.4) # La Nina
    enso_ma_ElN = np.ma.masked_less(enso_rand,0.4) # El Nino
# mask SST time series according to the ENSO index
    sst_ma_LaN = SSTa_d[~enso_ma_LaN.mask,:,:].copy()
    sst_ma_ElN = SSTa_d[~enso_ma_ElN.mask,:,:].copy()

    for ln in range(0,Y):
        for lt in range(0,X):
                    
            Pdays_LaN['n_p90_LaN'][rt,lt,ln] = np.count_nonzero(sst_ma_LaN[:,lt,ln] \
                                                                > sst_p90[lt,ln], \
                                                                axis=0)
            Pdays_LaN['n_p10_LaN'][rt,lt,ln] = np.count_nonzero(sst_ma_LaN[:,lt,ln] \
                                                                < sst_p10[lt,ln], \
                                                                axis=0)
                  
            Pdays_ElN['n_p90_ElN'][rt,lt,ln] = np.count_nonzero(sst_ma_ElN[:,lt,ln] \
                                                                > sst_p90[lt,ln], \
                                                                axis=0)
            Pdays_ElN['n_p10_ElN'][rt,lt,ln] = np.count_nonzero(sst_ma_ElN[:,lt,ln] \
                                                                < sst_p10[lt,ln], \
                                                                axis=0)
            
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



