'''

    compute the correlation maps of each of the 3 least oscillating 
    RCs from the SSA compared to the original time series at each locations

    Author: Eva C.
    Created: Oct 2018
    Last MOdif:

'''
######################
# libriries
######################
import numpy as np
import xarray as xr
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays
from scipy import signal
import pandas as pd
import time as time
import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

########################
# load data
########################
# EOF on the yearly mean (filtered) and 1x1 deg resolution grid
# and not weighted
fname_pcf = '../../ana/PotPred/vSZ/SSA/All_SSA_RCs_YearlyOnMonthly_NoWeight.nc'
tim_vec_f = xr.open_dataset(fname_pcf)['time']

# original time series at 1x1 deg resolution
fname_ssta = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/SSTa_monthly_trend.nc'
ssta = xr.open_dataset(fname_ssta)['SSTa']
lat = xr.open_dataset(fname_ssta)['lat']
lon = xr.open_dataset(fname_ssta)['lon']

# want to use the 3 top RCs from 6 PCs
maps = np.nan*np.zeros((6*3,len(lat),len(lon)))
pvalue = np.nan*np.zeros((6*3,len(lat),len(lon)))
mm=1 # mode
rc=1 # RC
t_key = time.time()
for ii in range(0,6*3):
    if rc<4:
        test = xr.open_dataset(fname_pcf)['RC_allPCs'][rc-1,:,mm-1].squeeze()
        corr, p = np.apply_along_axis(st.spearmanr,0,ssta,test)
        maps[ii,:,:] = corr
        pvalue[ii,:,:] = p
        rc = rc+1
    elif rc>3:
        mm = mm+1
        rc = 1
        test = xr.open_dataset(fname_pcf)['RC_allPCs'][rc-1,:,mm-1].squeeze()
        corr, p = np.apply_along_axis(st.spearmanr,0,ssta,test)
        maps[ii,:,:] = corr
        pvalue[ii,:,:] = p
        rc = rc+1

elapsed_key = time.time() - t_key
print('elapsed time to calc all correlation maps:', elapsed_key)

##########################################
# plotting
##########################################

lev = np.hstack((np.arange(-1,-0.1+0.1,0.1), \
                 np.arange(0.1,1+0.1,0.1)))
lev_cb = np.arange(-1,1+0.1,0.2)


plt.figure(figsize=(12,17))
plt.clf()
panel = 1
for ii in range(0,3*3):
    ax = plt.subplot(3,3,panel,projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines('50m', linewidth=0.8)
    ax.gridlines()
    plt.contourf(lon, lat, maps[ii,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
    cb=plt.colorbar(ticks=lev_cb,shrink=0.8)
    plt.contour(lon, lat, maps[ii,:,:], levels=[0], color='k')

    panel = panel+1

plt.figure(figsize=(12,17))
plt.clf()
panel = 1
for ii in range(3*3,6*3):
    ax = plt.subplot(3,3,panel,projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines('50m', linewidth=0.8)
    ax.gridlines()
    plt.contourf(lon, lat, maps[ii,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
    cb=plt.colorbar(ticks=lev_cb,shrink=0.8)
    plt.contour(lon, lat, maps[ii,:,:], levels=[0], color='k')

    panel = panel+1

#plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)
plt.show(block=False)


#### may want to use the p value to plot!!! seems to have a weak significance 
## which I think could make sense


