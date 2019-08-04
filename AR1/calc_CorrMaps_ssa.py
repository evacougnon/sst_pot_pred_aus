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
fname_pcf = '../tmp_data2plot/All_SSA_RCs_daily_NoWeight_365.nc'
tim_vec_f = xr.open_dataset(fname_pcf)['time']

# original time series at 1x1 deg resolution
fname_ssta = '../tmp_data2plot/SSTa_daily_trend_Aus_1deg.nc'
ssta = xr.open_dataset(fname_ssta)['SSTa']
#ssta = np.ma.masked_where(np.isnan(ssta)==True, ssta, copy=False)
# !!! when using the masked array, python wants the 2 input in the
# !!! st.spearmanr to be the same size, can't manage to get 
# !!! np/apply_along_axis to work! but we can assume that all the nans
# !!! at time step 0 are land ponts only
lat = xr.open_dataset(fname_ssta)['lat']
lon = xr.open_dataset(fname_ssta)['lon']

figfile1 = '../tmp_data2plot/corr_maps_1degAus_trend_PC1-4_RC1-3_daily_NoWeight_365_.eps'
#figfile2 = '/v_Munk_Drive/ecougnon/ana/PotPred/vSZ/SSA/corr_maps_1degAus_trend_PC4-56_RC1-3_daily_NoWeight_365.eps'

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
        maps[ii,:,:] = np.where(np.isnan(ssta[0,:,:])==True,np.nan,maps[ii,:,:])
        pvalue[ii,:,:] = p
        pvalue[ii,:,:] = np.where(np.isnan(ssta[0,:,:])==True,np.nan,pvalue[ii,:,:])
        rc = rc+1
    elif rc>3:
        mm = mm+1
        rc = 1
        test = xr.open_dataset(fname_pcf)['RC_allPCs'][rc-1,:,mm-1].squeeze()
        corr, p = np.apply_along_axis(st.spearmanr,0,ssta,test)
        maps[ii,:,:] = corr
        maps[ii,:,:] = np.where(np.isnan(ssta[0,:,:])==True,np.nan,maps[ii,:,:])
        pvalue[ii,:,:] = p
        pvalue[ii,:,:] = np.where(np.isnan(ssta[0,:,:])==True,np.nan,pvalue[ii,:,:])
        rc = rc+1

elapsed_key = time.time() - t_key
print('elapsed time to calc all correlation maps:', elapsed_key)

############################################################
# create masked array for significance (using the p value)
##############################################################
mask_sig = np.ones(pvalue.shape)
#for ii in range(0,6*3):
mask_sig = np.ma.masked_where(abs(maps)<pvalue,mask_sig)   
mask_sig = np.ma.masked_where(np.isnan(maps)==True,mask_sig)

##########################################
# plotting
##########################################
lev = np.hstack((np.arange(-0.8,-0.1+0.1,0.1), \
                 np.arange(0.1,0.8+0.1,0.1)))
lev_cb = np.arange(-0.8,0.8+0.1,0.2)

mpl.rcParams['hatch.linewidth'] = 0.5
mpl.rcParams['hatch.color'] = 'grey'

#plt.figure(figsize=(12,17))
my_dpi = 300
plt.figure(figsize=(4000/my_dpi,3500/my_dpi))
plt.clf()
panel=1
for ii in range(0,1*3):
    ax = plt.subplot(4,3,panel,projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines('50m', linewidth=0.8)
    ax.gridlines()
    plt.contourf(lon, lat, -maps[ii,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
    cb=plt.colorbar(ticks=lev_cb,shrink=0.8)
    plt.contour(lon, lat, -maps[ii,:,:], levels=[0], color='k')
    plt.contourf(lon, lat, mask_sig[ii,:,:],color='none', hatches = ['///'], \
                 edgecolr=0.6,alpha=0)
    ax.set_xlim([90, 180])
    ax.set_ylim([-55, 10])
    ax.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
    
    panel = panel+1

for ii in range(3,4*3):
    ax = plt.subplot(4,3,panel,projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines('50m', linewidth=0.8)
    ax.gridlines()
    plt.contourf(lon, lat, maps[ii,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
    cb=plt.colorbar(ticks=lev_cb,shrink=0.8)
    plt.contour(lon, lat, maps[ii,:,:], levels=[0], color='k')
    plt.contourf(lon, lat, mask_sig[ii,:,:],color='none', hatches = ['///'], \
                 edgecolr=0.6,alpha=0)
    ax.set_xlim([90, 180])
    ax.set_ylim([-55, 10])
    ax.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())

    panel = panel+1

plt.savefig(figfile1, bbox_inches='tight', format='eps', dpi=300)
'''
#plt.figure(figsize=(12,17))
plt.figure(figsize=(4000/my_dpi,3500/my_dpi))
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
    plt.contourf(lon, lat, mask_sig[ii,:,:],color='none', hatches = ['///'], \
                 edgecolr=0.6,alpha=0)
    ax.set_xlim([90, 180])
    ax.set_ylim([-55, 10])
    ax.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())

    panel = panel+1

plt.savefig(figfile2, bbox_inches='tight', format='eps', dpi=300)
plt.show(block=False)
'''



