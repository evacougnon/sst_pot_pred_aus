'''
    HadISST read and check mean of 1961-1990 with 
    1982-2005 
'''
# load required modules

import numpy as np
import xarray as xr
import pandas as pd

from scipy import io # load matlab file
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm
import cmocean

figfile ='/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SSTmean_hadisst_19611990_19822005.png'

# define indexes for lat lon of one pixel
lat_px_min = -46 # deg N
lat_px_max = -26 
lon_px_min = 135 #150 # deg E
lon_px_max = 174  

ds_hadisst = xr.open_dataset('/home/data/sst/HadISST_sst.nc')['sst']. \
                sel(latitude=slice(lat_px_max,lat_px_min), \
                    longitude=slice(lon_px_min,lon_px_max))

## do area-averaged (however, on a regular 1deg grid), depends if the 
# wanted area-averaged needs to be based per degree of m
# for 1961-1990
ds_hadisst_19611990 = ds_hadisst.sel(time=slice('1961-01-01','1990-12-31')). \
                                 mean(dim=('time'))
ds_hadisst_19822005 = ds_hadisst.sel(time=slice('1982-01-01','2005-12-31')). \
                                 mean(dim=('time'))

ds_hadisst_offset = ds_hadisst_19822005 - ds_hadisst_19611990

## plotting
domain = [lat_px_min, lon_px_min, lat_px_max, lon_px_max] 
domain_draw = [lat_px_min, lon_px_min, lat_px_max, lon_px_max]
dlat = 5 #30 #10
dlon = 10 #90 #30
llon, llat = np.meshgrid(ds_hadisst.longitude, ds_hadisst.latitude)
bg_col = '0.6'
cont_col = '1.0'
lev_temp = np.arange(10,26+2,2)
lev_diff = np.arange(-0.6,0.6+0.1,0.1)


plt.figure(figsize=(16,5))
plt.clf()
ax1 = plt.subplot(1,3,1)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, ds_hadisst_19611990, levels = lev_temp, \
             cmap=cmocean.cm.thermal)
cb=plt.colorbar(ticks=lev_temp, shrink=0.8)
plt.title('Mean SST 1961-1990', fontsize=16, y=1.04)

ax2 = plt.subplot(1,3,2)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, ds_hadisst_19822005, levels = lev_temp, \
             cmap=cmocean.cm.thermal)
cb=plt.colorbar(ticks=lev_temp, shrink=0.8)
plt.title('Mean SST 1982-2005', fontsize=16, y=1.04)

ax3 = plt.subplot(1,3,3)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, ds_hadisst_offset, levels = lev_diff, \
             cmap=cmocean.cm.balance)
cb=plt.colorbar(ticks=lev_diff, shrink=0.8)
plt.title('Difference', fontsize=16, y=1.04)

plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

#plt.show()









