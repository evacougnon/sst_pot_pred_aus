# plt std

import numpy as np
import xarray as xr

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm


# Load data

#fname = '/home/ecougnon/ana/SST_map_mhw_Aus'
fname = '/home/ecougnon/ana/SSTa_monthly_extremes_filter_Aus_20032017.nc'
fname_ = '/home/ecougnon/ana/SSTa_monthly_extremes_Aus.nc'
#SSTa_daily_filter_Aus_20032017.nc'
#SSTa_monthly_extremes_Aus.nc' #SSTa_daily_Aus.nc'
figfile = '/home/ecougnon/ana/SSTa_monthly_stats_filter_notfilter_std.png'

lat_map = xr.open_dataset(fname)['lat']
lon_map = xr.open_dataset(fname)['lon']
#SSTa  = xr.open_dataset(fname)['SSTa']
TMM = xr.open_dataset(fname)['TMM']
Tp90 = xr.open_dataset(fname)['Tp90']
Tp10 = xr.open_dataset(fname)['Tp10']

#1/4 degree resolution (not smoothed)
lat_ = xr.open_dataset(fname_)['lat']
lon_ = xr.open_dataset(fname_)['lon']
TMM_ = xr.open_dataset(fname_)['TMM']. \
          sel(time=slice('2003-01-01','2017-10-31'))
Tp90_ = xr.open_dataset(fname_)['Tp90']. \
          sel(time=slice('2003-01-01','2017-10-31'))
Tp10_ = xr.open_dataset(fname_)['Tp10']. \
          sel(time=slice('2003-01-01','2017-10-31'))


#TMM_std = np.nanstd(SSTa,axis=0)
TMM_std = np.nanstd(TMM,axis=0)
Tp90_std = np.nanstd(Tp90,axis=0)
Tp10_std = np.nanstd(Tp10,axis=0)

#1/4 degree resolution (not smoothed
TMM_std_ = np.nanstd(TMM_,axis=0)
Tp90_std_ = np.nanstd(Tp90_,axis=0)
Tp10_std_ = np.nanstd(Tp10_,axis=0)

# plot setting
domain = [-55, 90, 10, 180] #[-80, -180, 85, 180] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon_map, lat_map)
llon_, llat_ = np.meshgrid(lon_,lat_)
#llat, llon = np.meshgrid(lat_map, lon_map)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure(figsize=(15,5)) #(19,7))
plt.clf()

plt.subplot(1,2,1, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_, llat_)
#plt.contourf(lonproj, latproj, TMM_std, levels=np.arange(0,1.2+0.1,0.1), \
#             cmap=plt.cm.YlOrBr)
plt.contourf(lonproj, latproj, Tp90_std_ - Tp10_std_, \
             levels=np.arange(-0.2,0.2+0.05,0.05), \
             cmap=plt.cm.bwr)
#cb=plt.colorbar(ticks=np.arange(0,1.2+0.1,0.2),shrink=0.8)
cb=plt.colorbar(ticks=np.arange(-0.2,0.2+0.1,0.1),shrink=0.8)
cb.ax.tick_params(labelsize=14)
#plt.title('STD monthly SSTa -- based on 2003-2017', \
#          fontsize=14, y=1.02)
plt.title('STD(SSTa Mp90) - STD(SSTa Mp10): 2003-2017; 1/4deg', \
          fontsize=14, y=1.02)

plt.subplot(1,2,2, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, Tp90_std - Tp10_std, \
             levels=np.arange(-0.2,0.2+0.05,0.05), \
#(0,1.2+0.1,0.1), \
             cmap=plt.cm.bwr) #YlOrBr)
#cb=plt.colorbar(ticks=np.arange(0,1.2+0.1,0.2),shrink=0.5)
cb=plt.colorbar(ticks=np.arange(-0.2,0.2+0.1,0.1),shrink=0.8)
cb.ax.tick_params(labelsize=14)
plt.title('STD(SSTa Mp90) - STD(SSTa Mp10): 2003-2017; 1deg smoothed', \
          fontsize=14, y=1.02)


'''
plt.subplot(1,3,3, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, Tp90_std, levels=np.arange(0,1.2+0.1,0.1), \
             cmap=plt.cm.YlOrBr)
cb=plt.colorbar(ticks=np.arange(0,1.2+0.1,0.2),shrink=0.5)
cb.ax.tick_params(labelsize=14)
plt.title('STD monthly p90 SSTa', \
          fontsize=16, y=1.02)
'''

plt.savefig(figfile,bbox_inches='tight', format='png', dpi=300)
plt.show()

