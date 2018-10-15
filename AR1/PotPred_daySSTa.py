'''

    Compute potential predictability on time series
    of each point following von Storch and Zwiers's book
    or the Zheng,Frederiksen,Lou method
   
    same than PotPred_monthly but adapted for daily
    Dealing with leap year -- remove 31 December to keep
     the continuity of each chunk
 
    run time (on the monthly mean, but should be the same with daily as the 
    number of chunk matters, not the length of each chunk):
              about 45 minutes for the region around Oz using
                   35 chunks
              likely to be just under 2 hours for 140 chunks
                   (3 months chunk length)    
              under 1.5 hour for the 6 months chunk length
'''

# import libraries
import numpy as np
from scipy import signal
from scipy import io
from scipy import stats
import pandas as pd
import xarray as xr
import time as time

import sys
sys.path.insert(0,'../libraries/')
import eric_oliver as eo

# load useful information
# useful numbers
# tau -- number of months per chunk
tau = 365 #12 
# days -- lengths of the chunk
outfile = '/home/ecougnon/ana/PotPred/PP_SSTa_daily_1yr_vSZ_Aus.nc'
fname = '../../PotPred/SSTa_daily_Aus.nc'
#SSTa_monthly_extremes_Aus.nc'
#OTEs_NumDays_month_Aus_19822016.nc'
#fname = '/home/ecougnon/data/HadISST_sst.nc'
##############
# be sure they are anomalies, if not set check = 1
check = 0
###############
deg_res = 0.25 #1 #5
lat_min = -55 #-80 #-20 #-55
lat_max = 10 #80 #20 #10
lat_bin = 1 #4 * deg_res
lon_min = 90 #1 #160 #90
lon_max = 180 #360 #270 #180
lon_bin = 1 #4 * deg_res
lat = xr.open_dataset(fname)['lat']
#lat = lat.sel(lat=slice(lat_min,lat_max,lat_bin))
lon = xr.open_dataset(fname)['lon']
#lon = lon.sel(lon=slice(lon_min,lon_max,lon_bin))
#tim = xr.open_dataset(fname)['time']
#tim = tim.sel(time=slice('1871-01-01','2016-01-01'))
tim = pd.date_range('1982-01-01','2016-12-31')
# remove the last day of the year when leap year
# Need a function!!!....
tim = tim[tim !='1984-12-31']
tim = tim[tim !='1988-12-31']
tim = tim[tim !='1992-12-31']
tim = tim[tim !='1996-12-31']
tim = tim[tim !='2000-12-31']
tim = tim[tim !='2004-12-31']
tim = tim[tim !='2008-12-31']
tim = tim[tim !='2012-12-31']
tim = tim[tim !='2016-12-31']
SST = xr.open_dataset(fname)['SSTa']
SST = SST.sel(time=tim, lat=lat, lon=lon)
T_keys = ['SSTa']
# time period
MinYear = 1982 #1870 #1982
MaxYear = 2016 # 2016
NumYears = MaxYear - MinYear+1
NumMonths = NumYears*tau
# define the start/end indexes for each chunk
str_id = range(0,NumMonths,tau)
end_id = range(tau-1,NumMonths+1,tau)
NumChunk = len(str_id)

# allocate memory
#PP_keys = ['Var_interC','Var_noise','Var_slow','p']
PP_keys = ['Var_interC', 'Var_noise', 'p', 'F90', 'F95']
'''
describe the keys ..
Var_interC -- total inter chunk variance
Var_noise -- variance of the noise over all chunks
Var_slow -- variance ofthe potentially predictable component
p -- Potential inter-chunk predictability
'''
X = len(lon)
Y = len(lat)
var = len(PP_keys)
PP = xr.Dataset({'TMM':(('PP_keys','lat','lon'),np.zeros((var, Y, X)))}, \
                {'PP_keys': PP_keys, 'lat': lat, 'lon': lon})

t_key = time.time()

# the data should already be anomalies (detrended and deseasonned)
SST = np.array(SST)

if check==1: # do anomalies (deseasonned and detrended)
    # detrend the time series
    SST[SST<-100]=np.nan
    SST_flat = SST.reshape(len(tim),X*Y)
    dsst = np.empty((len(tim),X*Y))
    dsst.fill(np.nan)
    tt=np.arange(0,len(tim))
    for i in range(0,len(SST_flat[0,:])):
        valid = ~np.isnan(SST_flat[:,i])
        if (valid.any()==True):
            y = SST_flat[:,i]
            mean, trend, alpha = eo.trend(tt,y)
            dsst[:,i] = y - (tt*trend) -mean
#            dsst[:,i] = signal.detrend(SST_flat[valid,i], axis=0, type='linear')
        elif (valid.all()==False):
            dsst[:,i] = np.nan
        tmp_dsea, sea, beta = eo.deseason_harmonic(dsst[:,i],4,12)
        dsst[:,i] = np.squeeze(tmp_dsea[:,0])
    DSST = dsst.reshape(len(tim),Y,X)
    SST = DSST.copy()

## apply PotPred_ZhengFred function on a single time series
#PP[key][:,:,:] = np.apply_along_axis(eac.PotPred_ZhengFred,0, SST[:,:,:],tau, \
#                                     NumChunk, str_id,end_id)

# apply PotPred_vStorchZwiers function on a single time series
PP[key][:,:,:] = np.apply_along_axis(eac.PotPred_vStorchZwiers,0, \
                                     SST[:,:,:],tau, \
                                     NumChunk, str_id,end_id)

elapsed_key = time.time() - t_key
print('elapsed time for each key:', elapsed_key)

# saving data
PP.to_netcdf(outfile)
print('data saved')


    
'''
icnt_id = 0
for i in range(0,len(SST['lon'])):
    t_lon = time.time()
    print(icnt_id+1, 'of', len(SST['lon']))
    
    jcnt_id = 0
    for j in range(0,len(SST['lat'])):

        
        for key in T_keys:
#            t_key = time.time()
# checking that the time series is not empty -- detrend does not
# handle nans... assume that if nans in tmm, same for all the others
            if ((np.isnan(np.min(SST[key][:,j,i])) == False) & \
                (np.isnan(np.max(SST[key][:,j,i])) == False)):
# detrend data -- using the linear least squares fit
                dsst = signal.detrend(SST[key][:,j,i],type='linear')
        
# apply PotPred_ZhengFred function on a single time series
                [var1, var2, var3, var4] = eac.PotPred_ZhengFred(dsst,tau, \
                                                                 NumChunk, \
                                                                 str_id,end_id)
                PP[key][0,jcnt_id,icnt_id] = var1
                PP[key][1,jcnt_id,icnt_id] = var2
                PP[key][2,jcnt_id,icnt_id] = var3
                PP[key][3,jcnt_id,icnt_id] = var4
            
#            elapsed_key = time.time() - t_key
#            print('elapsed time for each key:', elapsed_key)

        jcnt_id = jcnt_id + 1

    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt_id = icnt_id + 1
'''


#'''
## climatology
##                t_clim = time.time()
#                clim = SST[key][:,j,i].groupby('time.month').mean('time')
#                sst_a = SST[key][:,j,i].groupby('time.month') - clim
##                elapsed_clim = time.time() - t_clim
##                print('elapsed time for each clim:', elapsed_clim)
## detrend data -- using the linear least squares fit
##                t_det = time.time()
#                dsst = signal.detrend(sst_a,type='linear')
##                elapsed_det = time.time() - t_det
##                print('elapsed time for each detrend:', elapsed_det)
#'''

