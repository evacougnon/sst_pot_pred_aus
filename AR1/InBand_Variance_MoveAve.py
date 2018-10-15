'''
	In Band variance following Terry's suggestion
	and Didier M.'s paper (2016)

	Created: Jul 2018
	Author: Eva C.

	STILL WORKING ON IT!! -- only using a moving average for now!!!
'''
##################
# import libraries
###################
import numpy as np
import xarray as xr
#from scipy import signal
from eofs.xarray import Eof
import dask.array as da
import bottleneck as bn

from xarray.core import dask_array_ops

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean

# allowing parallel calculations
#from dask.distributed import Client
import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries')
import xrft # detrend is possible there
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


###########################
# I/O file names
#########################
#outfile = 
figfile = '/home/ecougnon/ana/PotPred/InBand/'

############################
# define region and load data
# specify a lat and lon bin if subsetting!
# TMM: monthly temperature 
# Tp90: monthly 90th percentile
# TMX: monthly max
# Tp10: monthly 10th percentile
# TMN: monthly min
#
# 'HadISST' to use for the HadISST file
# already monthly SSTs
############################
WHICH = 'TMM' #'HadISST' #'TMM'
lat_min = -56 
lat_max = 10 
lat_bin = 1 
lon_min = 90 
lon_max = 180 
lon_bin = 1 
if WHICH != 'HadISST':
    fname = '/home/ecougnon/ana/SSTa_monthly_extremes_Aus.nc'
    deg_res = 0.25 
    lat = xr.open_dataset(fname)['lat'].sel(lat=slice(lat_min,lat_max,lat_bin))
    lon = xr.open_dataset(fname)['lon'].sel(lon=slice(lon_min,lon_max,lon_bin))
    tim = xr.open_dataset(fname)['time']
    SSTa_TMm = xr.open_dataset(fname)[WHICH].sel(time=tim, lat=lat, lon=lon)
elif WHICH=='HadISST':
    fname = '/home/data/sst/HadISST_sst.nc'
    deg_res = 1 
    lat = xr.open_dataset(fname)['latitude']. \
             sel(latitude=slice(lat_max,lat_min,lat_bin))
    lon = xr.open_dataset(fname)['longitude']. \
             sel(longitude=slice(lon_min,lon_max,lon_bin))
    tim = xr.open_dataset(fname)['time'].sel(time=slice('1871-01-01','2017-12-31'))
    ds_hadisst = xr.open_dataset(fname, chunks={'longitude':10, \
                                 'latitude':11})['sst']. \
                     sel(latitude=lat, longitude=lon, time=tim)
# detrend
    ds_hadisst = _apply_detrend(ds_hadisst,[0])
    clim_19611990 = ds_hadisst.sel(time=slice('1961-01-01', '1990-11-01')). \
                               groupby('time.month').mean(dim=('time'))
    clim_19822016 = ds_hadisst.sel(time=slice('1982-01-01', '2016-11-01')). \
                               groupby('time.month').mean(dim=('time'))

    SSTa_TMm = ds_hadisst.chunk({'time':12}).groupby('time.month') - clim_19822016

    SSTa_TMm.to_netcdf('/home/ecougnon/ana/SSTa_hadisst_Aus.nc')

#    SSTa_TMm = xr.open_dataset('/home/ecougnon/ana/SSTa_hadisst_Aus.nc')


SSTa_TMm_mean = SSTa_TMm.mean(dim='time')
SSTa_TMm_var = SSTa_TMm.var(dim='time')

############################
# apply a low pass filter
# boxcar
############################
if WHICH != 'HadISST':
    SSTa_TMm_f1 = SSTa_TMm.rolling(time=12, center=True).mean()
    SSTa_TMm_f2 = SSTa_TMm.rolling(time=12*2, center=True).mean()
    SSTa_TMm_f5 = SSTa_TMm.rolling(time=12*5, center=True).mean()
    SSTa_TMm_f10 = SSTa_TMm.rolling(time=12*10, center=True).mean()

elif WHICH=='HadISST':

#SSTa_TMm_f1 = dask_array_ops.dask_rolling_wrapper(bn.move_mean, \
#                                                  SSTa_TMm.chunk({'time':12}),12, \
#                                                  min_count=None, axis=0)

    SSTa_TMm_f1 = SSTa_TMm.chunk({'time':12}).rolling(time=12, center=True).mean()
    SSTa_TMm_f2 = SSTa_TMm.chunk({'time':12}).rolling(time=12*2, center=True).mean()
    SSTa_TMm_f5 = SSTa_TMm.chunk({'time':12}).rolling(time=12*5, center=True).mean()
    SSTa_TMm_f10 = SSTa_TMm.chunk({'time':12}).rolling(time=12*10, center=True).mean()

################################
# mean and varaiance
################################
SSTa_TMm_f1_mean = SSTa_TMm_f1.mean(dim='time')
SSTa_TMm_f1_var = SSTa_TMm_f1.var(dim='time')

SSTa_TMm_f2_mean = SSTa_TMm_f2.mean(dim='time')
SSTa_TMm_f2_var = SSTa_TMm_f2.var(dim='time')

SSTa_TMm_f5_mean = SSTa_TMm_f5.mean(dim='time')
SSTa_TMm_f5_var = SSTa_TMm_f5.var(dim='time')

SSTa_TMm_f10_mean = SSTa_TMm_f10.mean(dim='time')
SSTa_TMm_f10_var = SSTa_TMm_f10.var(dim='time')

#######################################
# Pot pred ratio with filtered data
#######################################
potpred_TMm_f1 = SSTa_TMm_f1_var/SSTa_TMm_var
potpred_TMm_f2 = SSTa_TMm_f2_var/SSTa_TMm_var
potpred_TMm_f5 = SSTa_TMm_f5_var/SSTa_TMm_var
potpred_TMm_f10 = SSTa_TMm_f10_var/SSTa_TMm_var

potpred_TMm_f1f2 = SSTa_TMm_f2_var/SSTa_TMm_f1_var
potpred_TMm_f2f5 = SSTa_TMm_f5_var/SSTa_TMm_f2_var
potpred_TMm_f5f10 = SSTa_TMm_f10_var/SSTa_TMm_f5_var


#############################
# Signal to noise ratio
#############################
SNR_TMm_f1 = (SSTa_TMm_mean-SSTa_TMm_f1_mean)/ \
             np.sqrt(SSTa_TMm_var-SSTa_TMm_f1_var)
SNR_TMm_f2 = (SSTa_TMm_f1_mean-SSTa_TMm_f2_mean)/ \
             np.sqrt(SSTa_TMm_f1_var-SSTa_TMm_f2_var)
SNR_TMm_f5 = (SSTa_TMm_f2_mean-SSTa_TMm_f5_mean)/ \
             np.sqrt(SSTa_TMm_f2_var-SSTa_TMm_f5_var)
SNR_TMm_f10 = (SSTa_TMm_f5_mean-SSTa_TMm_f10_mean)/ \
             np.sqrt(SSTa_TMm_f5_var-SSTa_TMm_f10_var)

'''
#####################################
# EOF analysis
#####################################
sst_var = SSTa_TMm_f10.dropna(dim='time', how='all')
coslat = np.cos(np.deg2rad(sst_var.coords['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(sst_var, weights=wgts)
# Retrieve the leading EOFs and the leading PC time series
# EOFs and PCs are scaled to unit variance 
# (divided by the square-root of their eigenvalues)
eofs = solver.eofs(neofs=8)
pcs = solver.pcs(npcs=8, pcscaling=1)
# variance fraction
vf = solver.varianceFraction(neigs=8)
# north test
NT_err = solver.northTest(neigs=8, vfscaled=True)

sname_eof = '/home/ecougnon/ana/PotPred/InBand/EOFs_TMM_f10.nc'
sst_var = sst_var.rename('SSTa_filter10yr')
eof_param = xr.merge([eofs, pcs, vf, NT_err, sst_var])

eof_param.to_netcdf(sname_eof)

print(vf)

'''
#####################################
# ploting
#####################################
# variances maps
################
cmap1 = cmocean.cm.speed
lev = np.arange(0,1+0.1,0.1)
lev2 = np.arange(0,0.5+0.1,0.1)

plt.figure(figsize=(5,3))
'''
ax1 = plt.subplot(1,3,1,projection = ccrs.PlateCarree())
SSTa_TMm_var.plot(ax=ax1,levels=lev,cmap=cmap1)
ax1.coastlines()
plt.title('SSTa_TMm_var')
ax2 = plt.subplot(1,3,2,projection = ccrs.PlateCarree())
SSTa_TMm_f1_var.plot(ax=ax2,levels=lev2,cmap=cmap1)
ax2.coastlines()
plt.title('SSTa_TMm_f1_var')
ax3 = plt.subplot(1,3,3,projection = ccrs.PlateCarree())
'''
ax3 = plt.subplot(1,1,1,projection = ccrs.PlateCarree())
potpred_TMm_f1.plot(ax=ax3,levels=lev,cmap=cmap1)
ax3.coastlines()
plt.title('var 1-yr filtered/var tot')
#plt.savefig(figfile + 'var_NotFilt_1yrFilt_' + WHICH + '.png', \
#            bbox_inches='tight', format='png', dpi=300)
plt.savefig(figfile + 'PotPred_1yrFilt_' + WHICH + '.png', \
            bbox_inches='tight', format='png', dpi=300)

plt.figure(figsize=(5,3))
'''
ax1 = plt.subplot(1,3,1,projection = ccrs.PlateCarree())
SSTa_TMm_f1_var.plot(ax=ax1,levels=lev2,cmap=cmap1)
ax1.coastlines()
plt.title('SSTa_TMm_f1_var')
ax2 = plt.subplot(1,3,2,projection = ccrs.PlateCarree())
SSTa_TMm_f2_var.plot(ax=ax2,levels=lev2,cmap=cmap1)
ax2.coastlines()
plt.title('SSTa_TMm_f2_var')
ax3 = plt.subplot(1,3,3,projection = ccrs.PlateCarree())
'''
ax3 = plt.subplot(1,1,1,projection = ccrs.PlateCarree())
potpred_TMm_f2.plot(ax=ax3,levels=lev,cmap=cmap1)
ax3.coastlines()
plt.title('var 2-yr filtered/var tot') #1-yr filtered')
plt.savefig(figfile + 'PotPred_2yrFilt_' + WHICH + '.png', \
            bbox_inches='tight', format='png', dpi=300)

plt.figure(figsize=(5,3))
'''
ax1 = plt.subplot(1,3,1,projection = ccrs.PlateCarree())
SSTa_TMm_f2_var.plot(ax=ax1,levels=lev2,cmap=cmap1)
ax1.coastlines()
plt.title('SSTa_TMm_f2_var')
ax2 = plt.subplot(1,3,2,projection = ccrs.PlateCarree())
SSTa_TMm_f5_var.plot(ax=ax2,levels=lev2,cmap=cmap1)
ax2.coastlines()
plt.title('SSTa_TMm_f5_var')
ax3 = plt.subplot(1,3,3,projection = ccrs.PlateCarree())
'''
ax3 = plt.subplot(1,1,1,projection = ccrs.PlateCarree())
potpred_TMm_f5.plot(ax=ax3,levels=np.arange(0,0.5+0.5,0.05),cmap=cmap1)
ax3.coastlines()
plt.title('var 5-yr filtered/var tot') #2-yr filtered')
plt.savefig(figfile + 'PotPred_5yrFilt_' + WHICH + '.png', \
            bbox_inches='tight', format='png', dpi=300)

#plt.figure(figsize=(17,3))
plt.figure(figsize=(5,3))
'''
ax1 = plt.subplot(1,3,1,projection = ccrs.PlateCarree())
SSTa_TMm_f5_var.plot(ax=ax1,levels=lev2,cmap=cmap1)
ax1.coastlines()
plt.title('SSTa_TMm_f5_var')
ax2 = plt.subplot(1,3,2,projection = ccrs.PlateCarree())
SSTa_TMm_f10_var.plot(ax=ax2,levels=lev2,cmap=cmap1)
ax2.coastlines()
plt.title('SSTa_TMm_f10_var')
ax3 = plt.subplot(1,3,3,projection = ccrs.PlateCarree())
'''
ax3 = plt.subplot(1,1,1,projection = ccrs.PlateCarree())
potpred_TMm_f10.plot(ax=ax3,levels=np.arange(0,0.5+0.5,0.05),cmap=cmap1)
ax3.coastlines()
plt.title('var 10-yr filtered/var tot') #5-yr filtered')
plt.savefig(figfile + 'PotPred_10yrFilt_' + WHICH + '.png', \
            bbox_inches='tight', format='png', dpi=300)

'''
###########
# SNR maps
###########
plt.figure(figsize=(21,3))
ax1 = plt.subplot(1,4,1,projection = ccrs.PlateCarree())
SNR_TMm_f1.plot(ax=ax1,cmap=cmap1)
ax1.coastlines()
plt.title('SNR_TMm_f1')
ax2 = plt.subplot(1,4,2,projection = ccrs.PlateCarree())
SNR_TMm_f2.plot(ax=ax2,cmap=cmap1)
ax2.coastlines()
plt.title('SNR_TMm_f2')
ax3 = plt.subplot(1,4,3,projection = ccrs.PlateCarree())
SNR_TMm_f5.plot(ax=ax3,cmap=cmap1)
ax3.coastlines()
plt.title('SNR_TMm_f5')
ax4 = plt.subplot(1,4,4,projection = ccrs.PlateCarree())
SNR_TMm_f10.plot(ax=ax4,cmap=cmap1)
ax4.coastlines()
plt.title('SNR_TMm_f10')
plt.savefig(figfile + 'SNR_InBand_' + WHICH + '.png', \
            bbox_inches='tight', format='png', dpi=300)



######
# EOFs
######
# plot North test:
plt.figure()
plt.errorbar(np.arange(1,8+1),vf*100, \
             yerr=NT_err[:8]*100,fmt='.')
plt.title('Eigenvalues spectrum -- North et al. 1982')
plt.xlabel('EOF modes')
plt.ylabel('Eigenvalue (%)')
plt.grid()
plt.savefig(figfile + 'NorthTest_' + WHICH + '_f10.png', \
            bbox_inches='tight', format='png', dpi=300)

# plot the 3 first EOF (can be changed if the north test
# shows more than 3 leading modes) with their PCs
clevs = np.arange(-0.008,0.008+0.01,0.001)

plt.figure(figsize=(17,7))
ax1 = plt.subplot(2,3,1,projection = ccrs.PlateCarree())
fill = eofs[0].plot.contourf(ax=ax1, levels=clevs, cmap=plt.cm.RdBu_r,
                             add_colorbar=True, transform=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, facecolor='w', edgecolor='k')
#cb = plt.colorbar(fill, orientation='horizontal')
ax1.set_title('EOF1', fontsize=16)

ax2 = plt.subplot(2,3,4)
pcs[:, 0].plot(color='k', linewidth=2)
ax2 = plt.gca()
#ax2.axhline(0, color='k')
ax2.set_ylim(-3, 3)
ax2.set_xlabel('Years')
ax2.set_ylabel('Normalized Units')
ax2.set_title('PC1 Time Series', fontsize=16)
plt.grid()

ax3 = plt.subplot(2,3,2,projection = ccrs.PlateCarree())
fill = eofs[1].plot.contourf(ax=ax3, levels=clevs, cmap=plt.cm.RdBu_r,
                             add_colorbar=True, transform=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, facecolor='w', edgecolor='k')
ax3.set_title('EOF2', fontsize=16)

ax4 = plt.subplot(2,3,5)
pcs[:, 1].plot(color='k', linewidth=2)
ax4 = plt.gca()
ax4.set_ylim(-3, 3)
ax4.set_xlabel('Years')
ax4.set_ylabel('Normalized Units')
ax4.set_title('PC2 Time Series', fontsize=16)
plt.grid()

ax5 = plt.subplot(2,3,3,projection = ccrs.PlateCarree())
fill = eofs[2].plot.contourf(ax=ax5, levels=clevs, cmap=plt.cm.RdBu_r,
                             add_colorbar=True, transform=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, facecolor='w', edgecolor='k')
ax5.set_title('EOF3', fontsize=16)

ax6 = plt.subplot(2,3,6)
pcs[:, 3].plot(color='k', linewidth=2)
ax6 = plt.gca()
ax6.set_ylim(-3, 3)
ax6.set_xlabel('Years')
ax6.set_ylabel('Normalized Units')
ax6.set_title('PC3 Time Series', fontsize=16)
plt.grid()

plt.savefig(figfile + 'EOFs_PCs_' + WHICH + '_f10.png', \
            bbox_inches='tight', format='png', dpi=300)


#plt.show()

'''





