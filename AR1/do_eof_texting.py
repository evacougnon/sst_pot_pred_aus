#######################################3
# implementing my own script for eof
# analysis
#########################################

# import libraries
import numpy as np
from scipy import linalg
import xarray as xr
import pandas as pd
import time as time

import dask.array

import sys
sys.path.insert(0,'../libraries/')
import eac_useful as eac

t_start = time.time()

# useful information
outfile_mtx = '../../ana/PotPred/EOF/CovMtx_Aus_trend_daily_1deg'
outfile_ts = '../../ana/PotPred/EOF/MtxDropna_Aus_trend_daily_1deg.nc'
outfile = '../../ana/PotPred/EOF/eof_Aus_daily_trend_1deg_12modes_monthly.nc'
fname = '../../ana/SSTa_daily_trend_Aus.nc'
#SSTa_daily_Aus.nc'
#SSTa_monthly_extremes_Aus.nc'
#SSTa_trend_monthly_Aus_19822016.nc'
#SSTa_with_season_monthly_Aus.nc'
deg_res = 1 #1 #5 # resolution wanted to avoid out of mem
lat_min = -55 
lat_max = 10 
lat_bin = 4 * deg_res # factor of 4 as the data are 1/4 degree
lon_min = 90 
lon_max = 180 
lon_bin = 4 * deg_res
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

SST_ = xr.open_dataset(fname)['SSTa'].sel(lat =slice(lat_min,lat_max,lat_bin), \
                                          lon=slice(lon_min,lon_max,lon_bin), \
                                           time=tim)
SST_ = SST_.chunk({'lat':13,'lon':18,'time':73}) 
print('chunking done')
print(SST_)

SST_ = SST_.transpose('time','lat','lon')
print('transpose done')
print(SST_)

# to do the EOF analysis. 
# Do a spatially weighted anomaly covariance matrix of a field
# The SSTs are already anomalies, then weight before the computation of EOFs: 
# Square-root of cosine of latitude
coslat = np.cos(np.deg2rad(SST_.coords['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]
SST_ = SST_/wgts

# flatten the lat/lon dimension
SST_ = SST_.stack(loc=('lat','lon'))
print('stacking in location done')
print(SST_)

SST_ = SST_.chunk({'time':5})

# averageeach chunk (annual mean here)
SST_m = SST_.resample(time='1MS').mean('time')
SST_ = SST_.groupby('time.year').mean('time').rename({'year':'time'})


t_key = time.time()

# drop land points
SST_ = SST_.dropna(dim='loc')
SST_m = SST_m.dropna(dim='loc')

elapsed_dropna = time.time() - t_key
print('elapsed time to drop nans', elapsed_dropna)

t_key = time.time()

# do the covariance matrix
cov = np.cov(SST_,rowvar=0)

elapsed_cov = time.time() - t_key
print('elapsed time to do the cov matrix', elapsed_cov)

'''

def covariance_gufunc(x, y):
    return ((x - x.mean(axis=0, keepdims=True))
            * (y - y.mean(axis=0, keepdims=True))).mean(axis=0)
def vector_norm(x, dim, ord=None):
    return xr.apply_ufunc(np.linalg.norm, x,
                          input_core_dims=[[dim]],
                          kwargs={'ord': ord, 'axis': -1})

xmean = ts_mtx.mean(dim='time')
cov = np.nansum(np.dot((ts_mtx - xmean),(ts_mtx - xmean)), \
                axis=0)/(ts_mtx.shape[0]-1)

'''

#input_core_dims=[['time']], 
#test = xr.apply_ufunc(np.cov, SST_, dask='allowed', kwargs={'rowvar':0})

### make the cov matrix a xarray to save along SST_
###COV = xr.Dataset({'cov':(('loc_1','loc_2'), np.zeros((

print('saving step 1 -- the covariance matrix!')
# saving
#SST_.unstack('loc').to_netcdf(outfile_ts)
#np.savez(outfile_mtx, cov = cov)


print('eigenvalues and vectors')
t_key = time.time()
#############################
# eigenvalues and vectors
############################
'''

EIGVAL = {}
keys_eigval = ['eigvalP','delta_eigval']
for key in keys_eigval:
    EIGVAL[key] = np.zeros(len(cov))

control = 1

'''

## calc eigenvalues and vectors
eigval, eigvec = linalg.eigh(cov)
# percentage of explained variance
eigvalP = eigval/np.sum(eigval) * 100
# Uses North et al equation 24 to see if eigenvalues (lambda) are 
# significantly separated
delta_eigval = eigvalP * np.sqrt(2/len(eigvalP))

control = 2

t_check = time.time()
# PCs, only the 12 first modes
modes=12
PCs = np.nan*np.zeros((modes,SST_m.shape[0]))
# transpose the SST_ before looping through!!!!
trans = SST_m.transpose('loc','time')
for pc in range(1,modes+1):
    PCs[-pc,:] = np.dot(eigvec[:,-pc],trans)

elapsed_pcs = time.time() - t_check
print('elapsed time to do the pcs', elapsed_pcs)

control = 3

# EOFmaps, only the 12 first modes
t_check = time.time()
modes=12
eof_maps = np.nan*np.zeros((modes,SST_.shape[1])) 
for ii in range(0,SST_.shape[1]): 
    eof_maps[:,ii] = eigvec[ii,-modes:]
elapsed_eofmaps = time.time() - t_check
print('elapsed time to do the eof maps', elapsed_eofmaps)

#PCs = np.apply_along_axis(np.multiply,1,SST_m,eof_maps[ii,:])
## finishe with a matrix of (time, loc) dimensions


EOFs = xr.Dataset({'eigvalP':(('modes'), eigvalP[-modes:]), \
                   'delta_eigval':(('modes'), delta_eigval[-modes:]), \
                   'PCs':(('modes','time'), PCs), \
                   'eof_maps':(('modes','loc'), eof_maps)}, \
                  coords = {'time':SST_m.time, \
                            'loc': SST_.lon.coords['loc'], \
                            'modes': np.arange(12,0,-1)})

# saving
EOFs.unstack('loc').to_netcdf(outfile)

elapsed_eof = time.time() - t_key
print('elapsed time to do the eofs', elapsed_eof)


elapsed_whole = time.time() - t_start
print('elapsed time since start', elapsed_whole)

