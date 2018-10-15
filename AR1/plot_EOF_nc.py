'''
    plot EOF maps with PCs from EOF analysis saved in netcdf format

  Author: Eva C.
  Created: Oct 2018 (from the plot_EOF.py based on npz files)
  Last Modif:

'''

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

# Load data
fname = '../../ana/PotPred/EOF/eof_Aus_daily_trend_1deg_12modes_monthly.nc' 
#eof_Aus_daily_1deg_12modes.nc'
fname_ppr = '../../ana/PotPred/PP_SSTa_daily_1yr_vSZ_Aus.nc'

figfile ='/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/EOFmaps_PCsm_5-8modes_SSTatrend.png'
figfile2 ='/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/NorthTest_SSTatrend.png'
#figfile_psd ='/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/PCsm_psd_5-8modes_SSTatrend.png'
#figfile_psdall ='/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/PCsm_psdall_5-8modes_SSTatrend_zoom.png'

lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
time_ts = xr.open_dataset(fname)['time']
MinYear=1982
MaxYear=2016

# what to plot
hmode = 5 # give the highest mode of 4 to be plotted

# load ppr
lat_map = xr.open_dataset(fname_ppr)['lat']
lon_map = xr.open_dataset(fname_ppr)['lon']
# ppr from vSZ if var_tot/var_fast
PP_vSZ = xr.open_dataset(fname_ppr)['TMM'][2,:,:]
PP_F90 = xr.open_dataset(fname_ppr)['TMM'][3,:,:]
'''
# to get ppr:
PP = (xr.open_dataset(fname_ppr)['TMM'][0,:,:] - \
      xr.open_dataset(fname_ppr)['TMM'][1,:,:]) \
     / xr.open_dataset(fname_ppr)['TMM'][0,:,:]
'''
#################################################
# NINO34
#################################################
file_nino34 = np.genfromtxt('/media/ecougnon/DATA/NewSverdrup/ecougnon/data/index/enso/nino34.txt', \
                            skip_header=1, skip_footer = 6)
str_nino34 = np.nonzero((file_nino34[:,0]>(MinYear-1)) \
                        & (file_nino34[:,0]<(MinYear+1)))
nino34_monthly = np.empty(len(time_ts))
k=0
for yy in np.arange(str_nino34[0][0],len(file_nino34[:,0])):
    for mm in np.arange(1,12+1):
        nino34_monthly[k] = file_nino34[yy,mm]
        k = k + 1
nino34_monthly = signal.detrend(nino34_monthly)
## plot the reconstructed time series for each mode
var_comp = nino34_monthly*10 # factor 10 to be readable on the figure
time_month = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
time_yr = pd.date_range('1982-01-01','2016-12-31',name='time',freq='12M')


#########################################################
## plotting 
########################################################
lev = np.hstack((np.arange(-0.07,-0.01+0.01,0.02), \
                 np.arange(0.01,0.07+0.01,0.02)))

plt.figure(figsize=(10,17)) 
plt.clf()

ax1 = plt.subplot(4,2,1,projection=ccrs.PlateCarree()) 
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
plt.contourf(lon, lat, xr.open_dataset(fname)['eof_maps'][-hmode,:,:], \
             levels = lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
#plt.contour(lon_map, lat_map, PP_vSZ, levels=[PP_F90], colors='C2')
plt.title('EOF' + str(hmode) +', ' \
          + str(np.around(xr.open_dataset(fname)['eigvalP'][-hmode].values,decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax2 = plt.subplot(4,2,3,projection=ccrs.PlateCarree()) 
ax2.add_feature(cfeature.LAND)
ax2.coastlines('50m', linewidth=0.8)
ax2.gridlines()
plt.contourf(lon, lat, xr.open_dataset(fname)['eof_maps'][-(hmode+1),:,:], \
             levels = lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
#plt.contour(lon_map, lat_map, PP_vSZ, levels=[PP_F90], colors='C2')
plt.title('EOF' + str(hmode+1) +', ' \
          + str(np.around(xr.open_dataset(fname)['eigvalP'][-(hmode+1)].values,decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax3 = plt.subplot(4,2,5,projection=ccrs.PlateCarree()) 
ax3.add_feature(cfeature.LAND)
ax3.coastlines('50m', linewidth=0.8)
ax3.gridlines()
plt.contourf(lon, lat, xr.open_dataset(fname)['eof_maps'][-(hmode+2),:,:], \
             levels = lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
#plt.contour(lon_map, lat_map, PP_vSZ, levels=[PP_F90], colors='C2')
plt.title('EOF' + str(hmode+2) +', ' \
          + str(np.around(xr.open_dataset(fname)['eigvalP'][-(hmode+2)].values,decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax4 = plt.subplot(4,2,7,projection=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND)
ax4.coastlines('50m', linewidth=0.8)
ax4.gridlines()
plt.contourf(lon, lat, xr.open_dataset(fname)['eof_maps'][-(hmode+3),:,:], \
             levels = lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
#plt.contour(lon_map, lat_map, PP_vSZ, levels=[PP_F90], colors='C2')
plt.title('EOF' + str(hmode+3) +', ' \
          + str(np.around(xr.open_dataset(fname)['eigvalP'][-(hmode+3)].values,decimals=1)) + \
          '%', fontsize=16, y=1.04)

ax5 = plt.subplot(4,2,2)
#r1, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-hmode,:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-hmode,:])
#ax_pc.annotate('r = {:.2f}'.format(r1), xy=(.1,.95), xycoords=ax_pc.transAxes)
plt.title('PC' + str(hmode) + ' with NINO3.4 index')
plt.grid()

ax6 = plt.subplot(4,2,4)
#r2, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-(hmode+1),:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-(hmode+1),:])
#ax_pc.annotate('r = {:.2f}'.format(r2), xy=(.1,.95), xycoords=ax_pc.transAxes)
#plt.title('PC' + str(hmode+1) + ' with NINO3.4 index')
plt.grid()

ax7 = plt.subplot(4,2,6)
#r3, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-(hmode+2),:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-(hmode+2),:])
#ax_pc.annotate('r = {:.2f}'.format(r3), xy=(.1,.95), xycoords=ax_pc.transAxes)
#plt.title('PC' + str(hmode+2) + ' with NINO3.4 index')
plt.grid()

ax8 = plt.subplot(4,2,8)
#r4, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-(hmode+3),:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-(hmode+3),:])
#ax_pc.annotate('r = {:.2f}'.format(r3), xy=(.1,.95), xycoords=ax_pc.transAxes)
#plt.title('PC' + str(hmode+3) + ' with NINO3.4 index')
plt.grid()

plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)
plt.show(block=False)

##################################################
## eigenvalues spectrum -- North et al. 1982
###################################################
plt.figure()
plt.errorbar(np.arange(0,12,1),xr.open_dataset(fname)['eigvalP'][:], \
             yerr=xr.open_dataset(fname)['delta_eigval'][:],fmt='.')
plt.title('Eigenvalues spectrum -- North et al. 1982')
plt.xlabel('EOF modes (1-12)')
plt.ylabel('Eigenvalue (%)')
plt.grid()
plt.savefig(figfile2, bbox_inches='tight', dpi=300)
plt.show(block=False)


###########################################
## spectral analysis on the PCs
###########################################
fs=1
f_dblside, pxx_dblside = np.apply_along_axis(signal.periodogram,0, \
                                             xr.open_dataset(fname)['PCs']. \
                                                transpose('time','modes'),fs, \
                                             window ='hanning',detrend=False, \
                                             return_onesided=False, \
                                             scaling='density')
window = 3
M = xr.open_dataset(fname)['PCs'].shape[0]
f_dblside_ = np.nan*np.zeros((len(f_dblside[:,0]),M))
pxx_dblside_ = np.nan*np.zeros((len(f_dblside[:,0]),M))
pxx_dblside_smooth = np.nan*np.zeros((len(f_dblside[:,0]),M))
for ii in range(0,M):
    id_sort = np.argsort(f_dblside[:,ii],axis=0)
    f_dblside_[:,ii] = f_dblside[id_sort,ii]
    pxx_dblside_[:,ii] = pxx_dblside[id_sort,ii]*2
    pxx_dblside_[int(len(pxx_dblside_)/2),ii] = \
                pxx_dblside_[int(len(pxx_dblside_)/2),ii]/2
    pxx_dblside_smooth[int(window/2):-int(window/2),ii]=  \
                       eac.movingaverage(pxx_dblside_[:,ii],window)


plt.figure(figsize=(10,17))
plt.clf()

ax5 = plt.subplot(4,2,1)
#r1, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-hmode,:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-hmode,:])
#ax_pc.annotate('r = {:.2f}'.format(r1), xy=(.1,.95), xycoords=ax_pc.transAxes)
plt.title('PC' + str(hmode) + ' with NINO3.4 index')
plt.grid()

ax6 = plt.subplot(4,2,3)
#r2, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-(hmode+1),:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-(hmode+1),:])
#ax_pc.annotate('r = {:.2f}'.format(r2), xy=(.1,.95), xycoords=ax_pc.transAxes)
#plt.title('PC' + str(hmode+1) + ' with NINO3.4 index')
plt.grid()

ax7 = plt.subplot(4,2,5)
#r3, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-(hmode+2),:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-(hmode+2),:])
#ax_pc.annotate('r = {:.2f}'.format(r3), xy=(.1,.95), xycoords=ax_pc.transAxes)
#plt.title('PC' + str(hmode+2) + ' with NINO3.4 index')
plt.grid()

ax8 = plt.subplot(4,2,7)
#r4, _ = st.spearmanr(xr.open_dataset(fname)['PCs'][-(hmode+3),:],var_comp)
ax_pc = plt.gca()
plt.plot(time_month,var_comp)
plt.plot(time_month,xr.open_dataset(fname)['PCs'][-(hmode+3),:])
#ax_pc.annotate('r = {:.2f}'.format(r3), xy=(.1,.95), xycoords=ax_pc.transAxes)
#plt.title('PC' + str(hmode+3) + ' with NINO3.4 index')
plt.grid()

ax1 = plt.subplot(4,2,2)
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-hmode], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-hmode],'*-')
plt.grid()

ax2 = plt.subplot(4,2,4)
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+1)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+1)],'*-')
plt.grid()

ax3 = plt.subplot(4,2,6)
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+2)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+2)],'*-')
plt.grid()

ax4 = plt.subplot(4,2,8)
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+3)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+3)],'*-')
plt.grid()


plt.savefig(figfile_psd, bbox_inches='tight', format='png', dpi=300)
plt.show(block=False)


 
