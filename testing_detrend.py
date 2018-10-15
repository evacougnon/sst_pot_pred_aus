'''
    testing a better way to detrend a multi dim array/xarray

    August 2018
    Eva Cougnon 
'''


import numpy as np
import xarray as xr
from scipy import signal
import time as time
import dask.array as da

import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries')
import xrft

lat_min = -35
lat_max = -30
lat_bin = 1
lon_min = 145
lon_max = 150
lon_bin = 1

# what is the input file(s)
# use WHICH = 'single' for a single input file
# WHICH = 'multiF' for multi files -- crashes easily due to memory issue
# when passing to a numpy array to loop through each location!

WHICH = 'multiF'

if WHICH == 'single':
    fname = '/home/data/sst/HadISST_sst.nc'
    deg_res = 1
    lat = xr.open_dataset(fname)['latitude']. \
             sel(latitude=slice(lat_max,lat_min,lat_bin))
    lon = xr.open_dataset(fname)['longitude']. \
             sel(longitude=slice(lon_min,lon_max,lon_bin))
    tim = xr.open_dataset(fname)['time'].sel(time=slice('1871-01-01','2017-12-31'))
    ds = xr.open_dataset(fname)['sst']. \
            sel(latitude=slice(lat_max,lat_min,lat_bin), \
                longitude=slice(lon_min,lon_max,lon_bin), \
                time=tim)
elif WHICH == 'multiF':
    header = '/home/data/sst/noaa_oi_v2/new_avhrr/'
    ds = xr.open_mfdataset(header + 'sst.day.mean.201?.v2.nc')['sst']. \
            sel(lat=slice(lat_min,lat_max), \
                lon=slice(lon_min,lon_max), \
                time=slice('1982-01-01','2018-08-06'))


#################
# new function to test!
# checkout: https://github.com/rabernat/xrft
# function copied into the libraries folder :)
################

def _apply_detrend(da, axis_num):
    """Wrapper function for applying detrending"""
    '''
	copied from xrft.py -- not sure why but the functions
	in xrft.py starting with '_' cannot be called directly!
	So I just copied here directly -- look for the reason later!!!!
    '''
    if da.chunks:
        func = xrft.detrend_wrap(xrft.detrendn)
        da = xr.DataArray(func(da.data, axes=axis_num),
                          dims=da.dims, coords=da.coords)
    else:
        if da.ndim == 1:
           da = xr.DataArray(sps.detrend(da),
                             dims=da.dims, coords=da.coords)
        else:
            da = xrft.detrendn(da, axes=axis_num)
        # else:
        #     raise ValueError("Data should be dask array.")
    return da


t_check = time.time()

# with this function, the array needs to be defined with chunks (as a dask array)
if WHICH == 'single':
    ds_chunk = ds.chunk({'longitude':1, 'latitude':1})
    test_new = _apply_detrend(ds_chunk,[0])
elif WHICH == 'multiF':
# with multi file, no need to chunk it as they are automatically chunked
# that's the way mfdataset works with dask array
    test_new = _apply_detrend(ds,[1])

elapsed = time.time() - t_check
print('elapsed time for the new function:', elapsed)




#################
# testing 
# going through a numpy array
################
#masked array not working with signal.detrend!!! #$%&^(#&@$$&#
#ds_hadisst_tmp = np.ma.masked_invalid(ds_hadisst_,copy=True)

# np.apply_along_axis with signal.detrend not working with land
# point as NaNs or even masked!....  
# test_old = np.apply_along_axis(signal.detrend,0,ds_hadisst_tmp)

# annoying looping way!
############################################
if WHICH == 'single':
    X = len(lon)
    Y = len(lat)
elif WHICH == 'multiF':
    X = len(ds.lon)
    Y = len(ds.lat)

t_check = time.time()

test_old = np.nan*np.zeros((len(ds.time),Y,X)) 

for ii in range(0,Y):
    for jj in range(0,X):
        valid = ~np.isnan(ds[:,ii,jj])
        if (valid.any()==True):
            test_old[:,ii,jj] = signal.detrend(np.squeeze(ds[:,ii,jj]))
        elif (valid.all()==False):
            test_old[:,ii,jj] = np.nan

elapsed = time.time() - t_check
print('elapsed time for the old way with', X*Y,'points:', elapsed)


#############
# do the diff the check they do the same thing!!
# can't compute test_new.compute() with NaNs -- which kind
# of makes sense as we shouldn't have to transform it in an 
# np.array! the whole ideais to keep using dask and xarray!
##############

diff_func = np.nan*np.zeros((len(tim),Y,X))
for ii in range(0,Y):
    for jj in range(0,X):
        valid = ~np.isnan(test_old[:,ii,jj])
        if (valid.any()==True):
            diff_func[:,ii,jj] = test_new[:,ii,jj] - test_old[:,ii,jj]
        elif (valid.all()==False):
            diff_func[:,ii,jj] = np.nan

print('max of the diffence new - old detrend calc:', np.nanmax(diff_func))
print('min of the diffence new - old detrend calc:', np.nanmin(diff_func))



