'''
    generate and save data to then use in plot_mapX.py
    data correspond to: number of days above the 90th percentile
    or 10th per pixels around Australia as well as the threshold 
    value for each pixel 

        firstly using the whole time series, then using the years 
          corresponding to a specific climate mode (will be region
          dependant!!!!)


	based on gendata_map_ensoYr_SSTa_significance.py 
	developed for the obs

	time estimated with 2500 random test of just under 8 hours

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
sname1 = header + 'SSTa_stats_nino_map_Aus_sig'
sname2 = header + 'SSTa_stats_nina_map_Aus_sig'

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

# number of random test
# in the previous metadata report with 35years of obsevations
# we had rtest 1000 random test, giving and equivalent of
# analysing 35000 years. To make it similar I chose 
# 2500 here as we have 15 years of data, giving an equivalent
# of 37500 years
rtest = 2500 

# allocate memory
Y = len(xt_ocean)
X = len(yt_ocean)
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
'''
SST_LaN = xr.Dataset({'Tp90':(('time','lon','lat'),np.zeros((rtest,X,Y))), \
                      'Tp10':(('time','lon','lat'),np.zeros((rtest,X,Y)))}, \
                     {'lat': yt_ocean, 'lon': xt_ocean, 'time': np.arange(0,rtest)})
SST_ElN = xr.Dataset({'Tp90':(('time','lon','lat'),np.zeros((rtest,X,Y))), \
                      'Tp10':(('time','lon','lat'),np.zeros((rtest,X,Y)))}, \
                     {'lat': yt_ocean, 'lon': xt_ocean, 'time': np.arange(0,rtest)})

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

##############################################
# randomising the time series rtest times
#############################################
for rt in range(0,rtest):
    rt_time = time.time()
    print(rt+1, 'of', rtest)
# randomly choose a starting index to 'reshuffle' the enso index time series
# so itr starts at any point of the enso time series and loop back to the 
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
# find p90 and p10 threshold as well as max and min
    SST_LaN['Tp90'][rt,:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                                  sst_ma_LaN[:,:,:],90)
    SST_LaN['Tp10'][rt,:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                                  sst_ma_LaN[:,:,:],10)
    SST_ElN['Tp90'][rt,:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                                  sst_ma_ElN[:,:,:],90)
    SST_ElN['Tp10'][rt,:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                                  sst_ma_ElN[:,:,:],10)
    del sst_ma_LaN
    del sst_ma_ElN
    elapsed_rt = time.time() - rt_time
    print('elapsed time for each random test:', elapsed_rt)

######################################################
# basic stats
######################################################
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






