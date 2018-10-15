# load required modules

import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
import time as time

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import marineHeatWaves as mhw

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

#figfile = '/home/ecougnon/ana/MHW_TASwest_area19822017.png'
'''
# define indexes for lat lon of one pixel
lat_px_min = -42.5 #-44 #-45 #-32 #-45 # deg N
lat_px_max = -42 # -37 #-28 #-37
lon_px_min = 148 #144 #147 #112 #147 # deg E
lon_px_max = 148.5 #150 #146 #155 #115 #155 

#time vector for the mhw function!!!
MinYear = 1982
MaxYear =2017
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,03).toordinal()+1)
'''

t_read = time.time()
header = '/home/data/sst/noaa_oi_v2/avhrr/'
ds = xr.open_mfdataset(header + '????/avhrr-only-v2.*.nc')['sst']. \
     sel(lat=slice(-42.5,-42), lon=slice(148,148.5))
elapsed_read = time.time() - t_read
print('elapsed time for read all the sst files: ', elapsed_read)

t_read = time.time()
ds_timeserie = xr.concat(ds, 'time')
elapsed_concat = time.time() - t_read
print('elapsed time for concatanating all the sst files along time: ', \
      elapsed_read)

# interm save
ds_timeserie.to_netcdf('/home/ecougnon/out.nc')

'''

# mhw framework onto the daily averaged over the little region
# !! keepin mind not area average!!
mhws, clim_h = mhw.detect(dtime, sst_obs, climatologyPeriod=[1982,2016], \
                          coldSpells=False)

# plot
# plot SSTa or daily time series
plt.figure(figsize=(13,13))
ax = plt.subplot(211)
plt.plot(tim_vec,ssta_obs)
plt.title('SSTa TAS area time series -- 42-44S 144-146E')
plt.grid()

ax = plt.subplot(212)
plt.plot(tim_vec, clim_h['seas'],'0.3')
plt.plot(tim_vec, sst_obs,'k')
plt.plot(tim_vec, clim_h['thresh'],'b')
plt.legend(['climatology','SST daily','threshold'])
plt.title('SST west TAS area -- 42-44S 144-146E -- with mhw threshold')
plt.grid()

#plt.savefig(figfile, bbox_inches='tight', dpi=300)

plt.show()
'''

'''
some info to plot with a bit more color, for a zoom in region

ax1.set_xlim([dtime_enso[0], dtime_enso[-1]])
ax1.set_ylim([-2, 3])
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
ax1.fill_between(dtime_enso, enso, where=enso >= 0.75, facecolor='red', alpha=0.5, \
                 interpolate=True)
ax1.fill_between(dtime_enso, enso, where=enso <= -0.75, facecolor='blue', alpha=0.5, \
                 interpolate=True
'''


