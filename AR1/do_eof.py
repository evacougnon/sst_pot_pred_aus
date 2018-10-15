#######################################3
# implementing my own script for eof
# analysis
#########################################

# import libraries
import numpy as np
from scipy import linalg
import xarray as xr
import time as time
import pandas as pd
import time as time

import dask.array

import sys
sys.path.insert(0,'../libraries/')
import eac_useful as eac

# useful information
tau = 365 #12 # length of each chunk for the potential pred
outfile_mtx = '/data/home/evac/NESP/PotPred/EOF/CovMtx_Aus_daily_1deg'
outfile = '/data/home/evac/NESP/PotPred/EOF/eof_Aus_daily_1deg_Modes1-8'
fname = '../../PotPred/SSTa_daily_Aus.nc'
#SSTa_monthly_extremes_Aus.nc'
#SSTa_trend_monthly_Aus_19822016.nc'
#SSTa_with_season_monthly_Aus.nc'
deg_res = 1 #5 # resolution wanted to avoid out of mem
lat_min = -55 #-80 #-20 #-55
lat_max = 10 #80 #20 #10
lat_bin = 1#4 * deg_res # factor of 4 as the data are 1/4 degree
lon_min = 90 #1 #160 #90
lon_max = 180 #360 #270 #180
lon_bin = 1 #4 * deg_res
lat = xr.open_dataset(fname)['lat']
lat = lat.sel(lat=slice(lat_min,lat_max,lat_bin))
lon = xr.open_dataset(fname)['lon']
lon = lon.sel(lon=slice(lon_min,lon_max,lon_bin))
#time_ts = xr.open_dataset(fname)['time']
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
## time period
MinYear = 1982
MaxYear = 2016
NumYears = MaxYear - MinYear+1
NumMonths = NumYears*tau

SST_ = xr.open_dataset(fname)['SSTa'].sel(lat =slice(lat_min,lat_max,lat_bin), \
                                          lon=slice(lon_min,lon_max,lon_bin), \
                                           time=tim)
SST_ = SST_.chunk({'lat':10,'lon':10})
print('chunking done')
print(SST_)

SST_ = SST_.transpose('time','lat','lon')
print('transpose done')
print(SST_)

# averageeach chunk (annual mean here)
SST_ = SST_.groupby('time.year').mean('time').rename({'year':'time'})

# Create an EOF solver to do the EOF analysis. 
# Do a spatially weighted anomaly covariance matrix of a field
# The SSTs are already anomalies, then weight before the computation of EOFs: 
# Square-root of cosine of latitude
coslat = np.cos(np.deg2rad(SST_.coords['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]

SST_ = SST_/wgts

print('Now trying to change the xarray as an np.array')
t_key = time.time()
SST = np.array(SST_)
elapsed_key = time.time() - t_key
print('elapsed time to get np.array', elapsed_key)

# reshape to have all the locations in the same dimension, time in the other
SST_flat = SST.reshape((t,X*Y))
lat_flat = np.repeat(np.copy(lat),len(lon))
lon_flat = np.tile(np.copy(lon),len(lat))
print('flatten location done')

# remove land points
ts_mtx = SST_flat[:,~np.isnan(SST_flat).all(0)]
lat_NoLand = lat_flat[~np.isnan(SST_flat[0,:])]
lon_NoLand = lon_flat[~np.isnan(SST_flat[0,:])]
print('removing nans done')

print('about to start cov matrix')
t_key = time.time()
# covariance matrix (every degree)
cov_mtx = np.cov(ts_mtx, rowvar=0)

elapsed_key = time.time() - t_key
print('elapsed time to calc. the np.cov (without nans):', elapsed_key)


print('saving step 1 -- the covariance matrix!')
# saving
np.savez(outfile_mtx, cov_mtx_chunk=cov_mtx_chunk, cov_noise=cov_noise, \
         lat_NoLand=lat_NoLand, lon_NoLand=lon_NoLand)

print('eigenvalues and vectors')


#############################
# eigenvalues and vectors
############################
EIGVAL = {}
keys_eigval = ['eigvalP','delta_eigval']
for key in keys_eigval:
    EIGVAL[key] = np.zeros(len(cov_mtx))

## calc eigenvalues and vectors
eigval, eigvec = linalg.eigh(cov_mtx)
# percentage of explained variance
EIGVAL['eigvalP'] = eigval_tot/np.sum(eigval_tot) * 100
# Uses North et al equation 24 to see if eigenvalues (lambda) are 
# significantly separated
EIGVAL['delta_eigval'] = EIGVAL['eigvalP'] * np.sqrt(2/len(EIGVAL['eigvalP']))

# PCs
PCs = np.zeros((len(cov_mtx),len(tim)))
for pc in range(1,len(cov_mtx)+1):
    PCs[pc-1,:] = np.dot(eigvec[:,-pc],np.transpose(ts_mtx))

# EOFmaps, only the 12 first modes
modes=12
eof_maps = np.zeros((modes,len(SST_flat[0,:])))
k=0
for ii in range(0,len(SST_flat[0,:])):
    if (np.isnan(SST_flat[0,ii]).all(0) == False):
       eof_maps[:,ii]= eigvec[k,-modes:]
       k += 1
    elif (np.isnan(SST_flat[0,ii]).all(0) == True):
       eof_maps[:,ii]= np.nan()

eof_maps = eof_maps.reshape(modes,Y,X)

# saving
np.savez(outfile, eof_maps=eof_maps, EIGVAL=EIGVAL, PCs=PCs, lat=lat, lon=lon, \
         tim=tim)




