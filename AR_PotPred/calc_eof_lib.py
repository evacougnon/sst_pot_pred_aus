#####################################
#
#  Do EOF analysis using the library to check with the TMM
#
#
#######################################3

# import libraries
import numpy as np
from scipy import linalg
import xarray as xr
import time as time
from eofs.xarray import Eof

from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import matplotlib
#import mpl_toolkits.basemap as bm

import sys
sys.path.insert(0,'../libraries/')

# load data
fname = '../../ana/SSTa_daily_Aus.nc'
#SSTa_trend_monthly_Aus_19822016.nc
#OTEs_NumDays_month_Aus_19822016.nc'
#SSTa_monthly_extremes_global.nc'
#SSTa_daily_Aus.nc'

#figfile = '/home/ecougnon/Desktop/WorkingFigures/EOFcov_SSTa_trend_EOFlib_yr.png'

deg_res = 1 #5 
lat_min = -55 #-80 #-20 #-55
lat_max = 10 #80 #20 #10
lat_bin = 1 #4 * deg_res  
lon_min = 90 #1 #160 #90
lon_max = 180 #360 #270 #180
lon_bin = 1 # 4 * deg_res
lat = xr.open_dataset(fname)['lat']
lat = lat.sel(lat=slice(lat_min,lat_max,lat_bin))
lon = xr.open_dataset(fname)['lon']
lon = lon.sel(lon=slice(lon_min,lon_max,lon_bin))

SST = xr.open_dataset(fname)['SSTa'] #['Pdays_p10'] #['TMM'] #['SSTa'] #['TMM']
SST = SST.sel(lat =slice(lat_min,lat_max,lat_bin), \
              lon=slice(lon_min,lon_max,lon_bin)).groupby('time.year'). \
          mean('time').rename({'year':'time'})

SST = SST.transpose('time','lat','lon')

#SST_ = SST.stacked(loc=('lat','lon'))
#SST_ = SST_.transpose('loc','time')

# Create an EOF solver to do the EOF analysis. 
# Do a spatially weighted anomaly covariance matrix of a field
# The SSTs are already anomalies, then weight before the computation of EOFs: 
# Square-root of cosine of latitude
coslat = np.cos(np.deg2rad(SST.coords['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(SST, weights=wgts, center=True)
#solver = Eof(SST_, weights=wgts, center=True)

# The Eof solver does a 'temporal' covariance matrix
# http://ajdawson.github.io/eofs/userguide/method.html?highlight=spatial%20covariance%20matrix


# center = True if we want to remove the mean; =False if no need to remove the mean
'''
If *True*, the mean along the first axis of *dataset* (the
            time-mean) will be removed prior to analysis. If *False*,
            the mean along the first axis will not be removed. Defaults
            to *True* (mean is removed).
            The covariance interpretation relies on the input data being
            anomaly data with a time-mean of 0. Therefore this option
            should usually be set to *True*. Setting this option to
            *True* has the useful side effect of propagating missing
            values along the time dimension, ensuring that a solution
            can be found even if missing values occur in different
            locations at different times.
'''
lambdas=solver.eigenvalues()
vf = solver.varianceFraction()
Nerror = solver.northTest(vfscaled=True)
pcs = solver.pcs() #(time, mode)
eofs = solver.eofsAsCovariance()
'''
plt.figure()
plt.subplot(3,2,1)
pcs[:, 0].plot()#color='b', linewidth=2)
ax = plt.gca()
ax.axhline(0, color='k')
ax.set_xlabel('Year')
ax.set_ylabel('PC1 amplitude')
plt.grid()
plt.subplot(3,2,2)
pcs[:, 1].plot()
ax = plt.gca()
ax.axhline(0, color='k')
ax.set_xlabel('Year')
ax.set_ylabel('PC2 amplitude')
plt.grid()
plt.subplot(3,2,3)
pcs[:, 2].plot()
ax = plt.gca()
ax.axhline(0, color='k')
ax.set_xlabel('Year')
ax.set_ylabel('PC3 amplitude')
plt.grid()
plt.subplot(3,2,4)
pcs[:, 3].plot()
ax = plt.gca()
ax.axhline(0, color='k')
ax.set_xlabel('Year')
ax.set_ylabel('PC4 amplitude')
plt.grid()
plt.subplot(3,2,5)
pcs[:, 4].plot()
ax = plt.gca()
ax.axhline(0, color='k')
ax.set_xlabel('Year')
ax.set_ylabel('PC5 amplitude')
plt.grid()
plt.subplot(3,2,6)
pcs[:, 5].plot()
ax = plt.gca()
ax.axhline(0, color='k')
ax.set_xlabel('Year')
ax.set_ylabel('PC6 amplitude')
#ax.set_title('PCs Time Series', fontsize=16)
plt.grid()
#plt.legend(['PC1','PC2','PC3','PC4','PC5','PC6'])
plt.show(block=False)

plt.figure()
plt.errorbar(np.arange(0,10),vf[:10]*100, \
             yerr=Nerror[:10]*100,fmt='.')
plt.title('Eigenvalues spectrum -- North et al. 1982')
plt.xlabel('EOF modes')
plt.ylabel('Eigenvalue (%)')
plt.grid()
plt.show(block=False)

'''


## plotting 
'''
domain = [int(lat_min), int(lon_min), int(lat_max), int(lon_max)]
domain_draw = [int(lat_min), int(lon_min), int(lat_max), int(lon_max)]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'
'''
lev = np.hstack((np.arange(-0.8,-0.1+0.1,0.1), np.arange(0.1,0.8+0.1,0.1)))
lev_ = np.arange(-0.8,0.8+0.2,0.2)
#lev = np.hstack((np.arange(-0.4,-0.01+0.01,0.02), np.arange(0.01,0.4+0.01,0.02)))

# EOF map

plt.figure(figsize=(12,11)) #12,11)) #7,11)) #(12,6)) # (10,8)
plt.clf()
ax1 = plt.subplot(3,2,1,projection=ccrs.PlateCarree())
'''
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
'''
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
plt.contourf(lon, lat, eofs[0,:,:], \
             levels= lev, \
             cmap=plt.cm.seismic, transform=ccrs.PlateCarree())
cb=plt.colorbar(ticks=lev_,shrink=0.9)
#cb.ax1.tick_params(labelsize=14)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('EOF 1, TMM ' + str(round(vf[0].values*100,1)) +'%', \
          fontsize=16, y=1.04)

ax2 = plt.subplot(3,2,2, projection=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND)
ax2.coastlines('50m', linewidth=0.8)
ax2.gridlines()
plt.contourf(lon, lat, eofs[1,:,:], levels=lev,  \
             cmap=plt.cm.seismic, transform=ccrs.PlateCarree())
cb=plt.colorbar(ticks=lev_,shrink=0.9)
#cb.ax2.tick_params(labelsize=14)
ax2.set_xlim([90, 180])
ax2.set_ylim([-55, 10])
ax2.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax2.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('EOF 2, TMM ' + str(round(vf[1].values*100,1)) +'%', \
          fontsize=16, y=1.04)

ax3 = plt.subplot(3,2,3,projection=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND)
ax3.coastlines('50m', linewidth=0.8)
ax3.gridlines()
plt.contourf(lon, lat, eofs[2,:,:], levels=lev,  \
             cmap=plt.cm.seismic, transform=ccrs.PlateCarree())
cb=plt.colorbar(ticks=lev_,shrink=0.9)
#cb.ax3.tick_params(labelsize=14)
ax3.set_xlim([90, 180])
ax3.set_ylim([-55, 10])
ax3.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax3.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('EOF 3, TMM ' + str(round(vf[2].values*100,1)) +'%', \
          fontsize=16, y=1.04)

ax4 = plt.subplot(3,2,4,projection=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND)
ax4.coastlines('50m', linewidth=0.8)
ax4.gridlines()
plt.contourf(lon, lat, eofs[3,:,:], levels=lev, \
             cmap=plt.cm.seismic, transform=ccrs.PlateCarree())
cb=plt.colorbar(ticks=lev_,shrink=0.9)
#cb.ax4.tick_params(labelsize=14)
ax4.set_xlim([90, 180])
ax4.set_ylim([-55, 10])
ax4.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax4.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('EOF 4, TMM ' + str(round(vf[3].values*100,1)) +'%', \
          fontsize=16, y=1.04)

ax5 = plt.subplot(3,2,5,projection=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND)
ax5.coastlines('50m', linewidth=0.8)
ax5.gridlines()
plt.contourf(lon, lat, eofs[4,:,:], levels=lev,  \
             cmap=plt.cm.seismic, transform=ccrs.PlateCarree())
#cb.ax5.tick_params(labelsize=14)
ax5.set_xlim([90, 180])
ax5.set_ylim([-55, 10])
ax5.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax5.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('EOF 5, TMM ' + str(round(vf[4].values*100,1)) +'%', \
          fontsize=16, y=1.04)

ax6 = plt.subplot(3,2,6,projection=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND)
ax6.coastlines('50m', linewidth=0.8)
ax6.gridlines()
plt.contourf(lon, lat, eofs[5,:,:], levels=lev, \
             cmap=plt.cm.seismic, transform=ccrs.PlateCarree())
#cb.ax6.tick_params(labelsize=14)
ax6.set_xlim([90, 180])
ax6.set_ylim([-55, 10])
ax6.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax6.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('EOF 6, TMM ' + str(round(vf[5].values*100,1)) +'%', \
          fontsize=16, y=1.04)

#plt.savefig(figfile, bbox_inches='tight')
plt.show(block=False)






