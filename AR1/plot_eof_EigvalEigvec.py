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
'''
fname_eigve = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/eigvec_NoWeight.nc'
fname_eigva = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/eigval_NoWeight.nc'
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PCs_YearlyOnMonthly_trend_NoWeight.nc'

'''
# EOF on the unfiltered (daily) SSTa and (1/4)x(1/4) deg resolution grid
fname_eigve = '../cdo_tool/eigvect_daily_1deg_NoWeight.nc'
fname_eigva = '../cdo_tool/eigval_daily_1deg_NoWeight.nc'
fname_pc = '../cdo_tool/PCs_daily_trend_1deg_NoWeight.nc'
#'
figfile = '../tmp_data2plot/EOF_PC_1-4_daily_NoWeight.eps'
lat = xr.open_dataset(fname_eigve, decode_times=False)['lat']
lon = xr.open_dataset(fname_eigve, decode_times=False)['lon']
MinYear = 1982
MaxYear = 2016
#tim_vector = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
#'
tim_vector = pd.date_range('1982-01-01','2016-12-31',name='time',freq='D')
# remove the last day of the year when leap year
# Need a function!!!....
tim_vector = tim_vector[tim_vector !='1984-12-31']
tim_vector = tim_vector[tim_vector !='1988-12-31']
tim_vector = tim_vector[tim_vector !='1992-12-31']
tim_vector = tim_vector[tim_vector !='1996-12-31']
tim_vector = tim_vector[tim_vector !='2000-12-31']
tim_vector = tim_vector[tim_vector !='2004-12-31']
tim_vector = tim_vector[tim_vector !='2008-12-31']
tim_vector = tim_vector[tim_vector !='2012-12-31']
tim_vector = tim_vector[tim_vector !='2016-12-31']
###############################
#'
eof_eigval = xr.open_dataset(fname_eigva, decode_times=False)['SSTa'].squeeze()
eof_var = eof_eigval[:].values/np.sum(eof_eigval).values * 100
eof_modes = xr.open_dataset(fname_eigve, decode_times=False)['SSTa'] #\
#            * np.sqrt(eof_eigval[:12]) # not sure about this normalisation!
eof_pcs =  xr.open_dataset(fname_pc)['SSTa']
hmode=1 # highest mode out of 4 to plot (only plot 4 EOF+PC

##########################################
# normalise to get the PCs in oC unit
# and corresponding EOF maps
##########################################
# Normalizing each Principal Component (dividing by max amplitude of that PC)
max_amp_pcs = np.max(np.abs(eof_pcs),axis=0) # for each PC
eof_pcs_norm = eof_pcs/(np.tile(np.copy(max_amp_pcs),(len(eof_pcs[:,0]), 1)))
# applyng reverse operation to eigenvectors (C matrix)
tmp = np.tile(np.copy(max_amp_pcs),(len(eof_modes[0,0,:]),len(eof_modes[0,:,0]),1))
tmp = np.transpose(tmp,(2,1,0))
eof_modes_norm = eof_modes*tmp

############################################
# Apply North Test
# Uses North et al equation 24 to see if eigenvalues (lambda) are 
# significantly separated
############################################
delta_eigval = eof_var * np.sqrt(2/len(eof_var))


####################################################
# plot
####################################################
'''
lev = np.hstack((np.arange(-4,-0.5+0.5,0.5), \
                 np.arange(0.5,4+0.5,0.5)))
lev_cb = np.arange(-4,4+1,1)
####################################################
lev = np.hstack((np.arange(-0.07,-0.01+0.01,0.02), \
                 np.arange(0.01,0.07+0.01,0.02)))
lev_cb = np.arange(-0.07,0.07+0.01,0.02)
'''

lev = np.hstack((np.arange(-2.4,-0.2+0.2,0.2), \
                 np.arange(0.2,2.4+0.2,0.2)))
lev_cb = np.arange(-2,2+0.5,0.5)


plt.figure(figsize=(10,12))
plt.clf()

ax1 = plt.subplot(2,2,1,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
plt.contourf(lon, lat, -eof_modes_norm[hmode-1,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev_cb,shrink=0.5)
cb.ax.tick_params(labelsize=12)
plt.contour(lon, lat, eof_modes_norm[hmode-1,:,:], levels=[0], color='k')
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
ax1.tick_params(labelsize=12)
plt.title('EOF' + str(hmode) +', ' \
          + str(np.round(eof_var[hmode-1],decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax2 = plt.subplot(2,2,2,projection=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND)
ax2.coastlines('50m', linewidth=0.8)
ax2.gridlines()
plt.contourf(lon, lat, eof_modes_norm[hmode,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev_cb,shrink=0.5)
cb.ax.tick_params(labelsize=12)
plt.contour(lon, lat, eof_modes_norm[hmode,:,:], levels=[0], color='k')
ax2.set_xlim([90, 180])
ax2.set_ylim([-55, 10])
ax2.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
ax2.tick_params(labelsize=12)
plt.title('EOF' + str(hmode+1) +', ' \
          + str(np.round(eof_var[hmode],decimals=1)) + \
          '%', fontsize=16, y=1.04)
ax3 = plt.subplot(2,2,3,projection=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND)
ax3.coastlines('50m', linewidth=0.8)
ax3.gridlines()
plt.contourf(lon, lat, eof_modes_norm[hmode+1,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev_cb,shrink=0.5)
cb.ax.tick_params(labelsize=12)
plt.contour(lon, lat, eof_modes_norm[hmode+1,:,:], levels=[0], color='k')
ax3.set_xlim([90, 180])
ax3.set_ylim([-55, 10])
ax3.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
ax3.tick_params(labelsize=12)
plt.title('EOF' + str(hmode+2) +', ' \
          + str(np.round(eof_var[hmode+1],decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax4 = plt.subplot(2,2,4,projection=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND)
ax4.coastlines('50m', linewidth=0.8)
ax4.gridlines()
plt.contourf(lon, lat, eof_modes_norm[hmode+2,:,:], levels=lev, cmap=plt.cm.RdYlBu_r)
cb=plt.colorbar(ticks=lev_cb,shrink=0.5)
cb.ax.tick_params(labelsize=12)
plt.contour(lon, lat, eof_modes_norm[hmode+2,:,:], levels=[0], color='k')
ax4.set_xlim([90, 180])
ax4.set_ylim([-55, 10])
ax4.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax4.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
ax4.tick_params(labelsize=12)
plt.title('EOF' + str(hmode+3) +', ' \
          + str(np.round(eof_var[hmode+2],decimals=1)) + \
          '%', fontsize=16, y=1.04)

'''
ax5 = plt.subplot(4,2,2)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs_norm[:,hmode-1],'k')
plt.title('PC' + str(hmode) + '')
plt.grid()

x6 = plt.subplot(4,2,4)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs_norm[:,hmode],'k')
plt.title('PC' + str(hmode+1) + '')
plt.grid()

ax7 = plt.subplot(4,2,6)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs_norm[:,hmode+1],'k')
plt.title('PC' + str(hmode+2) + '')
plt.grid()

ax8 = plt.subplot(4,2,8)
ax_pc = plt.gca()
plt.plot(tim_vector,eof_pcs_norm[:,hmode+2],'k')
plt.title('PC' + str(hmode+3) + '')
plt.grid()
'''

plt.savefig(figfile, bbox_inches='tight', format='eps', dpi=300)
plt.show(block=False)

'''
##################################################
## eigenvalues spectrum -- North et al. 1982
###################################################
plt.figure()
plt.errorbar(np.arange(1,13,1),eof_var[:12], yerr=delta_eigval[:12],fmt='.')
plt.title('Eigenvalues spectrum -- North et al. 1982')
plt.xlabel('EOF modes (1-12)')
plt.ylabel('Eigenvalue (%)')
plt.grid()
#plt.savefig(figfile2, bbox_inches='tight', dpi=300)
plt.show(block=False)
'''

