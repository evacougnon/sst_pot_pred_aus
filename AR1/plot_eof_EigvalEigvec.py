'''

	plot the EOF modes along with the PCs using the output
	from CDO (eigenvalues and eigenvectors)

	Author: Eva C.
	Created:; Oct 2018
	Last Modif.:

'''

# library
# import libraries
import numpy as np
import xarray as xr
import pandas as pd
from scipy import signal
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

import sys
sys.path.insert(0,'../libraries/')
import eac_useful as eac

import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#####################################################
# load data
#####################################################
# EOF on the yearly mean and 1x1 deg resolution grid
fname_eigve = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/outfile2_spa.nc'
fname_eigva = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/outfile1_spa.nc'
fname_pc1 = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PC_yearlyONmonthly_00000.nc'

'''
# EOF on the unfiltered (daily) SSTa and (1/4)x(1/4) deg resolution grid
fname_eigve = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/outfile2_spa.nc'
fname_eigva = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/outfile1_spa.nc'
'''
figfile = '/v_Munk_Drive/ecougnon/ana/PotPred/vSZ/EOF_yearly_cdo/EOF_PC_5-8_yearly.png'
lat = xr.open_dataset(fname_eigve)['lat']
lon = xr.open_dataset(fname_eigve)['lon']
tim = xr.open_dataset(fname_eigva)['time']
MinYear = 1982
MaxYear = 2016
tim_vector = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
eof_var = xr.open_dataset(fname_eigva)['SSTa'].squeeze()
eof_var = eof_var[:].values/np.sum(eof_var).values * 100
eof_modes = xr.open_dataset(fname_eigve)['SSTa']
eof_pcs = np.nan*np.zeros((12,xr.open_dataset(fname_pc1)['SSTa'].shape[0]))

for ii in range(0,12):
    if ii < 10:
        eof_pcs[ii,:] = xr.open_dataset('/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PC_yearlyONmonthly_0000' + str(ii) + '.nc')['SSTa'].squeeze()
    elif ii >=10:
        eof_pcs[ii,:] = xr.open_dataset('/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PC_yearlyONmonthly_000' + str(ii) + '.nc')['SSTa'].squeeze()

hmode=5 # highest mode out of 4 to plot (only plot 4 EOF+PC

####################################################
# plot
####################################################
lev = np.hstack((np.arange(-4,-0.5+0.5,0.5), \
                 np.arange(0.5,4+0.5,0.5)))
lev_cb = np.arange(-4,4+1,1)


plt.figure(figsize=(12,17))
plt.clf()

ax1 = plt.subplot(4,2,1,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
plt.contourf(lon, lat, eof_modes[hmode-1,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev,shrink=0.8)
plt.contour(lon, lat, eof_modes[hmode-1,:,:], levels=[0], color='k')
plt.title('EOF' + str(hmode) +', ' \
          + str(np.round(eof_var[hmode-1],decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax2 = plt.subplot(4,2,3,projection=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND)
ax2.coastlines('50m', linewidth=0.8)
ax2.gridlines()
plt.contourf(lon, lat, eof_modes[hmode,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev,shrink=0.8)
plt.contour(lon, lat, eof_modes[hmode,:,:], levels=[0], color='k')
plt.title('EOF' + str(hmode+1) +', ' \
          + str(np.round(eof_var[hmode],decimals=1)) + \
          '%', fontsize=16, y=1.04)
ax3 = plt.subplot(4,2,5,projection=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND)
ax3.coastlines('50m', linewidth=0.8)
ax3.gridlines()
plt.contourf(lon, lat, eof_modes[hmode+1,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev,shrink=0.8)
plt.contour(lon, lat, eof_modes[hmode+1,:,:], levels=[0], color='k')
plt.title('EOF' + str(hmode+2) +', ' \
          + str(np.round(eof_var[hmode+1],decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax4 = plt.subplot(4,2,7,projection=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND)
ax4.coastlines('50m', linewidth=0.8)
ax4.gridlines()
plt.contourf(lon, lat, eof_modes[hmode+2,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev,shrink=0.8)
plt.contour(lon, lat, eof_modes[hmode+2,:,:], levels=[0], color='k')
plt.title('EOF' + str(hmode+3) +', ' \
          + str(np.round(eof_var[hmode+2],decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax5 = plt.subplot(4,2,2)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs[hmode-1,:],'k')
plt.title('PC' + str(hmode) + '')
plt.grid()

x6 = plt.subplot(4,2,4)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs[hmode,:],'k')
plt.title('PC' + str(hmode+1) + '')
plt.grid()

ax7 = plt.subplot(4,2,6)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs[hmode+1,:],'k')
plt.title('PC' + str(hmode+2) + '')
plt.grid()

ax8 = plt.subplot(4,2,8)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs[hmode+2,:],'k')
plt.title('PC' + str(hmode+3) + '')
plt.grid()

plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)
#plt.show(block=False)


