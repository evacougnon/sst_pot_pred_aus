# load SSTs daily time series and estimate a percentage of the data
# above/below p90/p10 values during El Nino and La Nina yeats



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
ENSO = 'neg' # consider only the negative phases

header = '/home/ecougnon/ana/'
if ENSO == 'neg':
    outfile = header + 'NumDays_LaN_p90p10_Aus_filter.nc'
elif ENSO == 'pos':
    outfile = header + 'NumDays_ElN_p90p10_Aus_filter.nc'

################################
# define the region
################################
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
# number of days for 10% of the total number of days
tim_10th = len(dtime)/10
# Load data
# Daily SSTa
SSTa_d = xr.open_dataset(header + 'SSTa_daily_filter_Aus_20032017.nc')['SSTa']
lat_map = SSTa_d.lat
lon_map = SSTa_d.lon
# threshold values (from the daily SSTa time series)
sst_p10 = np.apply_along_axis(np.nanpercentile,0,SSTa_d,10)
sst_p90 = np.apply_along_axis(np.nanpercentile,0,SSTa_d,90)

#################################
# allocate memory
#################################
Y = len(lat_map)
X = len(lon_map)
'''
n_p90_LaN -- number of days above p90 corresponding to La Nina year
n_p10_LaN -- number of days below p10 corresponding to La Nina year
p_p90_LaN -- percentage of number of days above p90 corresponding to La Nina year
p_p10_LaN -- percentage of number of days below p10 corresponding to La Nina year
'''
if ENSO == 'neg':
    Pdays = xr.Dataset({'n_p90_LaN':(('lat','lon'),np.zeros((Y, X))), \
                        'n_p10_LaN':(('lat','lon'),np.zeros((Y, X))), \
                        'p_p90_LaN':(('lat','lon'),np.zeros((Y, X))), \
                        'p_p10_LaN':(('lat','lon'),np.zeros((Y, X)))}, \
                       {'lat': lat_map, 'lon': lon_map})
elif ENSO == 'pos':
    Pdays = xr.Dataset({'n_p90_ElN':(('lat','lon'),np.zeros((Y, X))), \
                        'n_p10_ElN':(('lat','lon'),np.zeros((Y, X))), \
                        'p_p90_ElN':(('lat','lon'),np.zeros((Y, X))), \
                        'p_p10_ElN':(('lat','lon'),np.zeros((Y, X)))}, \
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
dssta_oisst = np.array(ssta_oisst)
#dssta_oisst = np.apply_along_axis(signal.detrend,0,ssta_oisst, type='linear')
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

if ENSO == 'pos':
# SSTa for positive ENSO
    sst_ElN = np.array(SSTa_d)
    sst_ElN[nino34_n_id,:,:] = np.nan
    sst_ElN[nino34_neu_id,:,:] = np.nan
elif ENSO == 'neg':
# SSTa for negative ENSO
    sst_LaN = np.array(SSTa_d)
    sst_LaN[nino34_p_id,:,:] = np.nan
    sst_LaN[nino34_neu_id,:,:] = np.nan

for ln in range(0,X):
    for lt in range(0,Y):
        if ENSO == 'neg':
            Pdays['n_p90_LaN'][lt,ln] = np.count_nonzero(sst_LaN[:,lt,ln] \
                                                         > sst_p90[lt,ln], \
                                                         axis=0)
            Pdays['n_p10_LaN'][lt,ln] = np.count_nonzero(sst_LaN[:,lt,ln] \
                                                         < sst_p10[lt,ln], \
                                                         axis=0)
            Pdays['p_p90_LaN'][lt,ln] = Pdays['n_p90_LaN'][lt,ln]/tim_10th
            Pdays['p_p10_LaN'][lt,ln] = Pdays['n_p10_LaN'][lt,ln]/tim_10th 
        
        elif ENSO == 'pos':
            Pdays['n_p90_ElN'][lt,ln] = np.count_nonzero(sst_ElN[:,lt,ln] \
                                                         > sst_p90[lt,ln], \
                                                         axis=0)
            Pdays['n_p10_ElN'][lt,ln] = np.count_nonzero(sst_ElN[:,lt,ln] \
                                                         < sst_p10[lt,ln], \
                                                         axis=0)
            Pdays['p_p90_ElN'][lt,ln] = Pdays['n_p90_ElN'][lt,ln]/tim_10th
            Pdays['p_p10_ElN'][lt,ln] = Pdays['n_p10_ElN'][lt,ln]/tim_10th
        

## save files
Pdays.to_netcdf(outfile)


###############################################################
# significance test
#############################################################
# saving file
sname = header + 'NumDays_ElN_p90p10_Aus_95sign_2500_filter'
sname2 = header + 'NumDays_LaN_p90p10_Aus_95sign_2500_filter'
# number of random test
rtest = 2500
####################################
# allocate memory
####################################
#n_p90_LaN -- number of days above p90 corresponding to La Nina year
#n_p10_LaN -- number of days below p10 corresponding to La Nina year


Pdays_LaN = xr.Dataset({'n_p90_LaN':(('time','lat','lon'),np.zeros((rtest,Y,X))), \
                        'n_p10_LaN':(('time','lat','lon'),np.zeros((rtest,Y,X)))}, \
                       {'lat': lat_map, 'lon': lon_map, 'time':np.arange(0,rtest)})

Pdays_ElN = xr.Dataset({'n_p90_ElN':(('time','lat','lon'),np.zeros((rtest,Y,X))), \
                        'n_p10_ElN':(('time','lat','lon'),np.zeros((rtest,Y,X)))}, \
                       {'lat': lat_map, 'lon': lon_map, 'time':np.arange(0,rtest)})


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

    for ln in range(0,X):
        for lt in range(0,Y):

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


