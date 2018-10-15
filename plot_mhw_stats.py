# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import time as time
from datetime import date

from scipy import io # load matlab file
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm
import cmocean

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo
import marineHeatWaves as mhw
import eac_useful as eac

figfile ='/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SST_TasSeaBox_MHWtimeseriers_grid.png'

#####################
# define some conditions
##########################
# whether of not turing on the area-averaged (area-weighted)
# or just a spatial average 9no weight
what = 'A' # 'S'

# define indexes for lat lon of one pixel
lat_px_min = -46 #-45 #-46 # deg N
lat_px_max = -26 #-37 #-26
lon_px_min = 135 #139 #147 #135 #150 # deg E
lon_px_max = 174 #147 #155 #174

####################################################
## define the offset by using the HadISST dataset
####################################################
ds_hadisst = xr.open_dataset('/home/data/sst/HadISST_sst.nc')['sst']. \
                sel(latitude=slice(lat_px_max,lat_px_min), \
                    longitude=slice(lon_px_min,lon_px_max), \
                    time=slice('1871-01-01','2017-12-31'))
# do spatial or area-averaged of the SST
if what=='S':
    ds_hadisst_mean = ds_hadisst.mean(dim=('longitude','latitude'))
elif what=='A':
    ds_hadisst_all = xr.open_dataset('/home/data/sst/HadISST_sst.nc')['sst']. \
                        sel(latitude=slice(lat_px_max+1,lat_px_min-1), \
                             longitude=slice(lon_px_min-1,lon_px_max+1), \
                             time=slice('1871-01-01','2017-12-31'))
    lon = ds_hadisst_all.longitude
    lat = ds_hadisst_all.latitude
# get the distance between each grid points in meters    
    dx, dy = eac.dxdy(lat,lon)
    
    mask_land = np.ones((len(ds_hadisst.latitude),len(ds_hadisst.longitude)))
    tmp = np.squeeze(ds_hadisst[0,:,:])
    mask_land = np.ma.masked_where(np.isnan(tmp)==True, mask_land)
    area = dx * dy * mask_land
    weight = area/np.nansum(area)
    tot_weight = weight.sum()
    ds_hadisst_mean = ((ds_hadisst*weight).sum(dim=('longitude','latitude')))


ds_hadisst_clim = ds_hadisst_mean.sel(time=slice('1961-01-01','1990-12-31'))
tim_vec_hadisst = pd.date_range('1871-01-01','2017-12-31',name='time',freq='M')
## do area-averaged (however, on a regular 1deg grid), depends if the 
# wanted area-averaged needs to be based per degree of m
# for 1961-1990
ds_hadisst_19611990 = ds_hadisst.sel(time=slice('1961-01-01','1990-12-31')). \
                                 mean(dim=('time','longitude','latitude'))
ds_hadisst_19822005 = ds_hadisst.sel(time=slice('1982-01-01','2005-12-31')). \
                                 mean(dim=('time','longitude','latitude'))
ds_hadisst_offset = ds_hadisst_19822005 - ds_hadisst_19611990

#################################################
# load the correcponding SSTs
################################################3
header = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/'
if what == 'S':
    ds_oisst_d = xr.open_dataset(header + 'TasSea_Oliver2017_bdr.nc')['sst']
elif what == 'A':
    ds_oisst_d = xr.open_dataset(header + \
                                 'TasSea_Box_AreaAve_updated_May.nc')['sst']
#'TasSea_135_AreaAve_masked_updated.nc')['sst']

ds_oisst_m = ds_oisst_d.resample('1MS', dim='time', how='mean')


####################################
## apply the MHW code
####################################
#time vector for the mhw function!!!
# needs to be generated wth datetime format!
MinYear = 1982
MaxYear = 2018 # warninng, so far finishes in Feb!!
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
# warning finishes in Jan 2018 not Dec
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,5,19).toordinal()+1)
#dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,31).toordinal()+1)
dates = [date.fromordinal(tt.astype(int)) for tt in dtime]

sst_mean_d = np.array(ds_oisst_d)
dT = np.array(ds_hadisst_offset)
mhws, clim = mhw.detect(dtime, sst_mean_d[:len(dtime)], \
                        climatologyPeriod=[1982,2005], \
                        alternateClimatology=[dtime, sst_mean_d - dT])

######################################
# plotting from: 
# https://github.com/ecjoliver/marineHeatWaves/blob/master/docs/example_synthetic.ipynb
########################################
ev = np.argmax(mhws['intensity_max']) # Find largest event
print('Maximum intensity:', mhws['intensity_max'][ev], 'deg. C')
print('Average intensity:', mhws['intensity_mean'][ev], 'deg. C')
print('Cumulative intensity:', mhws['intensity_cumulative'][ev], 'deg. C-days')
print('Duration:', mhws['duration'][ev], 'days')
print('Start date:', mhws['date_start'][ev].strftime("%d %B %Y"))
print('End date:', mhws['date_end'][ev].strftime("%d %B %Y"))


plt.figure(figsize=(15,4))
ax=plt.subplot(111)
'''
plt.subplot(2,1,1)
# Plot SST, seasonal cycle, and threshold
plt.plot(dates, clim['seas'], '0.5')
plt.plot(dates, sst_mean_d, 'k')
plt.plot(dates, clim['thresh'], 'b')
plt.title('SST (black), seasonal climatology (gray), \
          threshold (blue), detected MHW events (shading)')
plt.xlim(['2015-01-01','2018-05-31']) #dates[0], dates[-1])
plt.ylim(sst_mean_d.min()-0.5, sst_mean_d.max()+0.5)
plt.ylabel(r'SST [$^\circ$C]')
plt.subplot(2,1,2)
'''
# Find indices for all ten MHWs before  n until the event of interest 
# wich is the most rcent event in this case, and shade accordingly
for ev0 in np.arange(ev-10, ev+1, 1):
    t1 = np.where(dtime==mhws['time_start'][ev0])[0][0]
    t2 = np.where(dtime==mhws['time_end'][ev0])[0][0]
    plt.fill_between(dates[t1:t2+1], sst_mean_d[t1:t2+1], \
                     clim['thresh'][t1:t2+1], \
                     color=(1,0.6,0.5))
# Find indices for MHW of interest and shade accordingly
t1 = np.where(dtime==mhws['time_start'][ev])[0][0]
t2 = np.where(dtime==mhws['time_end'][ev])[0][0]
plt.fill_between(dates[t1:t2+1], sst_mean_d[t1:t2+1], clim['thresh'][t1:t2+1], \
                 color='r')
# Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
plt.plot(dates, sst_mean_d, 'k', linewidth=2)
plt.plot(dates, clim['seas'], '0.5', linewidth=2)
plt.plot(dates, clim['thresh'], 'b', linewidth=2)
plt.legend(['SST','seasonal climatology','threshold'], \
           loc=3)
plt.xlim(['2015-01-01','2018-05-31']) 
ax.xaxis.set_major_locator(matplotlib.dates.YearLocator())
ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator())
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
#plt.xticks(np.arange(dates[-(36*3+140)],dates[-1],30))
#plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
plt.ylim(clim['seas'].min() - 1, \
         clim['seas'].max() + mhws['intensity_max'][ev]) # + 0.5)
plt.ylabel(r'SST [$^\circ$C]')
plt.grid()

plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

plt.show()


'''
plt.figure(figsize=(15,7))
# Duration
plt.subplot(2,2,1)
evMax = np.argmax(mhws['duration'])
plt.bar(range(mhws['n_events']), mhws['duration'], width=0.6, \
        color=(0.7,0.7,0.7))
#plt.bar(range(mhws['date_start'].strftime("%d %B %Y")), mhws['duration'], width=0.6, \
#        color=(0.7,0.7,0.7))

plt.bar(evMax, mhws['duration'][evMax], width=0.6, color=(1,0.5,0.5))
plt.bar(ev, mhws['duration'][ev], width=0.6, edgecolor=(1,0.,0.), \
        color='none')
plt.xlim(0, mhws['n_events'])
plt.ylabel('[days]')
plt.title('Duration')
# Maximum intensity
plt.subplot(2,2,2)
evMax = np.argmax(mhws['intensity_max'])
plt.bar(range(mhws['n_events']), mhws['intensity_max'], width=0.6, \
        color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['intensity_max'][evMax], width=0.6, color=(1,0.5,0.5))
plt.bar(ev, mhws['intensity_max'][ev], width=0.6, edgecolor=(1,0.,0.), \
        color='none')
plt.xlim(0, mhws['n_events'])
plt.ylabel(r'[$^\circ$C]')
plt.title('Maximum Intensity')
# Mean intensity
plt.subplot(2,2,4)
evMax = np.argmax(mhws['intensity_mean'])
plt.bar(range(mhws['n_events']), mhws['intensity_mean'], width=0.6, \
        color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['intensity_mean'][evMax], width=0.6, color=(1,0.5,0.5))
plt.bar(ev, mhws['intensity_mean'][ev], width=0.6, edgecolor=(1,0.,0.), \
        color='none')
plt.xlim(0, mhws['n_events'])
plt.title('Mean Intensity')
plt.ylabel(r'[$^\circ$C]')
plt.xlabel('MHW event number')
# Cumulative intensity
plt.subplot(2,2,3)
evMax = np.argmax(mhws['intensity_cumulative'])
plt.bar(range(mhws['n_events']), mhws['intensity_cumulative'], width=0.6, \
        color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['intensity_cumulative'][evMax], width=0.6, color=(1,0.5,0.5))
plt.bar(ev, mhws['intensity_cumulative'][ev], width=0.6, edgecolor=(1,0.,0.), \
        color='none')
plt.xlim(0, mhws['n_events'])
plt.title(r'Cumulative Intensity')
plt.ylabel(r'[$^\circ$C$\times$days]')
plt.xlabel('MHW event number')

#plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

plt.show()

'''

