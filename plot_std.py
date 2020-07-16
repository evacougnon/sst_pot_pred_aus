# plt std or variance

import numpy as np
import xarray as xr

from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Load data

fname = "/home/ecougnon/Documents/PotPred_paper/data/SSTa_daily_Aus_shift_time_dim.nc" #SST_daily_trend_Aus.nc" #SSTa_daily_Aus.nc"
figfile = "/home/ecougnon/Documents/PotPred_paper/figures/SSTa_daily_var_shift_time_file_05deg_color.eps"

lat_min = -55
lat_max = 10
lat_bin = 2
lon_min = 90
lon_max = 180
lon_bin = 2

lat_map = xr.open_dataset(fname)['lat'].sel(lat=slice(lat_min,lat_max,lat_bin))
lon_map = xr.open_dataset(fname)['lon'].sel(lon=slice(lon_min,lon_max,lon_bin))
tim = xr.open_dataset(fname)['time']
SSTa = xr.open_dataset(fname)['SSTa'].sel(time=tim, lat=lat_map, lon=lon_map)

TMM_std = SSTa.var("time")   #np.nanvar(SSTa,axis=0)

# plot setting
#domain = [-55, 90, 10, 180] #[-80, -180, 85, 180] #[-55, 90, 10, 180]
#domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
#dlat = 10 #30 #10
#dlon = 30 #90 #30
#llon, llat = np.meshgrid(lon_map, lat_map)
#llat, llon = np.meshgrid(lat_map, lon_map)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure(figsize=(19,7)) #15,5)) #(19,7))
plt.clf()
ax1=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#TMM_std.plot()
plt.contourf(lon_map, lat_map, TMM_std, levels=np.arange(0,1.6+0.1,0.1), \
             cmap=plt.cm.afmhot_r)
cb=plt.colorbar(ticks=np.arange(0,1.6+0.1,0.2),shrink=0.8)
cb.ax.tick_params(labelsize=14)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('variance daily SSTa', \
          fontsize=14, y=1.02)

plt.savefig(figfile,bbox_inches='tight', format='eps', dpi=300)
plt.show()

