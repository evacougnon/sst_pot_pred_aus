# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import datetime
from datetime import date
import time as time

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import marineHeatWaves as mhw

from scipy import io # load matlab file
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm
import cmocean
import matplotlib.patches as patches

#figfile = '/home/ecougnon/ana/MHW_TASeast_area_OliverEtAl2017.png'

# Running SSTa or the MHW script:
TEST = 'SSTa'

####
# 1982-2016
#####

# define indexes for lat lon of one pixel
lat_px_min = -46 #-44 #-44 #-45 #-32 #-45 # deg N
lat_px_max = -26 #-42 # -37 #-28 #-37
lon_px_min = 150 #144 #50 #145 #148 #144 #147 #112 #147 # deg E
lon_px_max = 174 #146 #50.5 #145 #148.5 #150 #146 #155 #115 #155 
# time vector 
tim_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='D')
tim_str = np.where(tim_vec == np.datetime64(datetime.datetime(1982,1,1)))
tim_end = np.where(tim_vec == np.datetime64(datetime.datetime(2016,12,31)))

header_obs = '/home/data/sst/noaa_oi_v2/avhrr/'
matobj_tmp = io.loadmat(header_obs + \
                        'timeseries/avhrr-only-v2.ts.0001.mat')
lat_tmp = matobj_tmp['lat']
lon_tmp = matobj_tmp['lon']
res = 0.25 # resolution of the data -- quarter degree
# find the closest index from the lat/lon defined above
lat_px_min_id = np.nonzero((lat_tmp>(lat_px_min-res)) & \
                           (lat_tmp<(lat_px_min)))
lat_px_max_id = np.nonzero((lat_tmp>(lat_px_max-res)) & \
                           (lat_tmp<(lat_px_max)))
lon_px_min_id = np.nonzero((lon_tmp>(lon_px_min-res)) & \
                           (lon_tmp<(lon_px_min)))
lon_px_max_id = np.nonzero((lon_tmp>(lon_px_max-res)) & \
                           (lon_tmp<(lon_px_max)))
# load time series from the given point
sst_obs = np.empty((len(np.arange(lon_px_min_id[0][0], \
                    lon_px_max_id[0][0]+1)), \
                    len(np.arange(lat_px_min_id[0][0], \
                    lat_px_max_id[0][0]+1)), len(tim_vec)))
sst_obs.fill(np.nan)
k=0
for i in range(lon_px_min_id[0][0],lon_px_max_id[0][0]+1):
    print(i)
    matobj = io.loadmat(header_obs + 'timeseries/avhrr-only-v2.ts.' \
                        + str(i).zfill(4) + '.mat')
    sst_obs[k,:,:] = matobj['sst_ts'][lat_px_min_id[0][0]:lat_px_max_id[0][0]+1, \
                                      tim_str[0][0]:tim_end[0][0]+1]
    k = k+1
sst_obs = np.nanmean(np.nanmean(sst_obs,0),0)

######
# 2017
#####
t_read = time.time()
header = '/home/data/sst/noaa_oi_v2/avhrr/'
ds_2017 = xr.open_mfdataset(header + '2017/avhrr-only-v2.2017*.nc')['sst']. \
           sel(lat=slice(lat_px_min,lat_px_max), lon=slice(lon_px_min,lon_px_max)). \
           squeeze('zlev')
elapsed_read = time.time() - t_read
print('elapsed time for read all the sst files: ', elapsed_read)

######
# 2018
######
t_read = time.time()
ds_2018 = xr.open_dataset(header + '2018/sst.day.mean.2018.v2.nc')['sst']. \
          sel(lat=slice(lat_px_min,lat_px_max), lon=slice(lon_px_min,lon_px_max))

t_read = time.time()
ds_timeserie = xr.concat(ds, 'time').squeeze()
elapsed_concat = time.time() - t_read
print('elapsed time for concatanating all the sst files along time: ', \
      elapsed_read)




if TEST =='SSTa'
# monthly mean


'''
# SSTa
fname = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
lat_obs = xr.open_dataset(fname)['lat']
lat_obs = lat_obs.sel(lat=slice(lat_px_min,lat_px_max))
lon_obs = xr.open_dataset(fname)['lon']
lon_obs = lon_obs.sel(lon=slice(lon_px_min,lon_px_max))
tim_obs = xr.open_dataset(fname)['time']
tim_obs = tim_obs.sel(time=slice('1982-01-01','2016-12-31'))
ssta_obs = xr.open_dataset(fname)['SSTa']
ssta_obs = ssta_obs.sel(time=tim_obs, lat=lat_obs,lon=lon_obs)
ssta_obs = np.nanmean(np.nanmean(ssta_obs,0),0)
'''

elif TEST=='MHW':
########
# MHW
#########
    ds_timeserie = ds_timeserie.mean(dim=('lon','lat'))
#time vector for the mhw function!!!
    MinYear = 1982
    MaxYear =2017
    NumYears = MaxYear-MinYear+1
    MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
    NumDays = 365*NumYears + MaxNumLeapYear
    dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,31).toordinal()+1)

    SST = np.empty((len(dtime)))
    SST.fill(np.nan)
    SST[0:len(sst_obs)] = sst_obs
    SST[len(sst_obs):len(sst_obs)+len(ds_timeserie)] = ds_timeserie

    t_read = time.time()
# mhw framework onto the daily averaged over the little region
# !! keepin mind not area average!! but the data points are on a regular 
# quarter degree grid! should be fine
    mhws, clim_h = mhw.detect(dtime, SST, climatologyPeriod=[1982,2005]) #, \
#                          coldSpells=False)
    elapsed_mhw = time.time() - t_read
    print('elapsed time for: ', elapsed_mhw)

#######
# plot
##########
# plot SSTa or daily time series
plt.figure(figsize=(13,13))
'''
ax = plt.subplot(211)
plt.plot(tim_vec,ssta_obs)
plt.title('SSTa TAS point time series -- 42S 148E')
plt.grid()
'''
time_vec_2017 =  pd.date_range('1982-01-01','2017-12-31',name='time',freq='D')
ax = plt.subplot(211)
plt.plot(time_vec_2017, clim_h['seas'],'0.3')
plt.plot(time_vec_2017, SST,'k')
plt.plot(time_vec_2017, clim_h['thresh'],'b')
plt.legend(['climatology','SST daily','threshold'])
ax.set_xlim(['1982-01-01','2017-12-03'])
ax.set_ylim([10, 22])
ax.fill_between(time_vec_2017, SST, where=SST>= clim_h['thresh'], \
                facecolor='red', alpha=0.5, \
                interpolate=True)
plt.title('SST E TAS  area (same than Oliver et al. 2017) -- 37-45S 147-155E -- with mhw threshold')
plt.grid()

ax = plt.subplot(212)
plt.plot(time_vec_2017, clim_h['seas'],'0.3')
plt.plot(time_vec_2017, SST,'k')
plt.plot(time_vec_2017, clim_h['thresh'],'b')
plt.legend(['climatology','SST daily','threshold'])
ax.set_xlim(['2016-01-01','2017-12-04'])
ax.set_ylim([10, 22])
ax.fill_between(time_vec_2017, SST, where=SST>= clim_h['thresh'], \
                facecolor='red', alpha=0.5, \
                interpolate=True)
plt.title('SST E TAS area (same than Oliver et al. 2017) -- 37-45S 147-155E -- with mhw threshold')
plt.grid()


#plt.savefig(figfile, bbox_inches='tight', dpi=300)

plt.show()


##########
# ploting the map from one specific day
# with the square of the western region choosen
##############
snapS = xr.open_dataset(header + '2017/avhrr-only-v2.20171121_preliminary.nc')['sst']. \
        sel(lat=slice(-47,-40), lon=slice(141,152)).squeeze()
lat_map = xr.open_dataset(header + '2017/avhrr-only-v2.20171121_preliminary.nc')['lat']. \
        sel(lat=slice(-47,-40)).squeeze()
lon_map = xr.open_dataset(header + '2017/avhrr-only-v2.20171121_preliminary.nc')['lon']. \
        sel(lon=slice(141,152)).squeeze()

# plot setting
domain = [-47, 141, -40, 152] #[-80, 0, 85, 360] #[-55, 90, 10, 180]
domain_draw = [-47, 141, -40, 152] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 1 #30 #10
dlon = 2 #90 #30
llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure(figsize=(10,8)) #12,6)) # (10,8)
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, snapS, levels=np.arange(12,19+1,1), \
             cmap=cmocean.cm.haline)
cb=plt.colorbar(ticks=np.arange(12,19+1,1),shrink=0.9)
cb.ax.tick_params(labelsize=14)
ax.add_patch(patches.Rectangle(-44,144,2,2,fill=False))
plt.title('SST from 21 Novembre 2017', \
          fontsize=16, y=1.08)
#plt.savefig(figfile,bbox_inches='tight', dpi=300)
plt.show()








