'''
RMSE based on the difference between member001 and the ensemble mean
'''

# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import time as time
from sklearn.metrics import mean_squared_error

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_useful as eac

# load data
header_out = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'

# define the region of interest
lat_min = -55
lat_max = 10
# !!! LONGITUDE in the model!!! from -280 to 80 deg E
lon_min = -270 #90 when deg east from 0 to 360
lon_max = -180 #180 when deg east from 0 to 360
yt_ocean = xr.open_dataset(gname)['yt_ocean']
yt_ocean = yt_ocean.sel(yt_ocean=slice(lat_min,lat_max))
xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
xt_ocean = xt_ocean.sel(xt_ocean=slice(lon_min,lon_max))

# do not use the year 2002 as it is spinning up
MinYear = 2003
MaxYear = 2017

# time vectors
tim_day_vec = pd.date_range('2003-01-01','2017-12-19',name='time',freq='D')


SST_ens_ = xr.open_mfdataset(header + 'ocean_daily_SST_enkf9_mem001_20??-20??.nc') \
                            ['temp'].sel(yt_ocean=yt_ocean, \
                                         xt_ocean=xt_ocean).squeeze('st_ocean')
SST_ens = np.array(SST_ens_)
SST_ens = SST_ens.reshape((len(tim_day_vec),len(yt_ocean)*len(xt_ocean)))

SST_mem_ = xr.open_mfdataset(header + 'ocean_daily_SST_enkf9_mem002_20??-20??.nc') \
                            ['temp'].sel(yt_ocean=yt_ocean, \
                                         xt_ocean=xt_ocean).squeeze('st_ocean')
SST_mem = np.array(SST_mem_)
SST_mem = SST_mem.reshape((len(tim_day_vec),len(yt_ocean)*len(xt_ocean)))

# assume that the nans are only for the land mask -- set them to zero
# it should't impact the rmse calc as they'll all be at the same value
SST_ens[np.isnan(SST_ens)==True] = 0
SST_mem[np.isnan(SST_mem)==True] = 0

rmse_ = np.empty(len(SST_ens[0,:]))
rmse_.fill(np.nan)
for ll in range(0,len(SST_ens[0,:])):
    rmse_[ll] = mean_squared_error(SST_ens[:,ll],SST_mem[:,ll])

rmse = rmse_.reshape((len(yt_ocean),len(xt_ocean)))

## plotting
domain = [-55, -270, 10, -180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
domain_draw = [-55, -270, 10, -180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon_obs, llat_obs = np.meshgrid(xt_ocean, yt_ocean)
llon_mdl, llat_mdl = np.meshgrid(xt_ocean, yt_ocean)
bg_col = '0.6'
cont_col = '1.0'

plt.figure(figsize=(11,11)) #12,11)) #7,11)) #(12,6)) # (10,8)
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_mdl, llat_mdl)
plt.contourf(lonproj, latproj, rmse, \
             levels= np.arange(0,1+0.1,0.1), \
             cmap=plt.cm.viridis)
cb=plt.colorbar(ticks=np.arange(0,1+0.1,0.1),shrink=0.9)
plt.title('RMSE ensemble-one member', fontsize=16, y=1.04)


plt.savefig(header_out + 'rmse_ensemble_mem002.png', \
            bbox_inches='tight', format='png', dpi=300)
plt.show()


