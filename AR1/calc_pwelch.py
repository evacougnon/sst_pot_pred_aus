'''
        Use the pwelch method to estimate an appropriate
        time scale to apply with the PotPred methods
 
        Created: Aug 2018
        Author: Eva C.

        STILL WORKING ON IT!!
'''
##################
# import libraries
###################
import numpy as np
import xarray as xr
from scipy import signal
#import dask.array as da
#import bottleneck as bn

#from xarray.core import dask_array_ops

import matplotlib.pyplot as plt

'''
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
'''

###########################
# I/O file names
#########################
#outfile = 
figfile = '/home/ecougnon/ana/PotPred/InBand/'

############################
# define region and load data
# specify a lat and lon bin if subsetting!
# TMM: monthly temperature 
# Tp90: monthly 90th percentile
# TMX: monthly max
# Tp10: monthly 10th percentile
# TMN: monthly min
#
# 'HadISST' to use for the HadISST file
# already monthly SSTs
############################
WHICH = 'TMM' #'HadISST' #'TMM'
lat_min = -30 #-56
lat_max = -20 #10
lat_bin = 1
lon_min = 105 #90
lon_max = 115 #180
lon_bin = 1
#if WHICH != 'HadISST':
#lat_wanted = -5
#lon_wanted = 170
fname = '/home/ecougnon/ana/SSTa_monthly_extremes_Aus.nc'
deg_res = 0.25
lat = xr.open_dataset(fname)['lat'].sel(lat=slice(lat_min,lat_max,lat_bin))
lon = xr.open_dataset(fname)['lon'].sel(lon=slice(lon_min,lon_max,lon_bin))
tim = xr.open_dataset(fname)['time']
SSTa_TMm = xr.open_dataset(fname)[WHICH].sel(time=tim, lat=lat, lon=lon)

#test = np.squeeze(SSTa_TMm[:,0,0])
#test =  xr.open_dataset(fname)[WHICH].sel(lat=lat_wanted, lon=lon_wanted, \
#                                          method='nearest')
test = SSTa_TMm.mean(dim=('lat','lon'))

#sampling rate -- mmonthy
fs = 1#2 #/((365/12)*24*3600) # cycle per day (1/(365/12))
# second per month: (365/12)*24*3600
f_, Pxx_ = signal.periodogram(test,fs, window='hanning', \
                              detrend=False, scaling='spectrum') #density')

f_0, Pxx_den_0 = signal.welch(test, fs, window='hanning', \
                              detrend=False, scaling='spectrum') #density')
# same than a normal periodogram 
f_1, Pxx_den_1 = signal.welch(test, fs, window='hanning', nperseg=len(test), \
                              detrend=False, scaling='spectrum') #density')
# trying different size of segments and overlap len(test)/2 = 2 segments
f_2, Pxx_den_2 = signal.welch(test, fs, window='hamming', nperseg=len(test)/2, \
                              detrend=False, scaling='spectrum') #density')
f_2no, Pxx_den_2no = signal.welch(test, fs, window='hamming', nperseg=len(test)/2, \
                              noverlap=0, detrend=False, scaling='spectrum')
f_4, Pxx_den_4 = signal.welch(test, fs, window='hamming', nperseg=len(test)/4, \
                              detrend=False, scaling='spectrum') #density')
f_4no, Pxx_den_4no = signal.welch(test, fs, window='hamming', nperseg=len(test)/4, \
                              noverlap=0, detrend=False, scaling='spectrum')

#################
# plotting
#################
plt.figure()
plt.plot(f_,Pxx_)
plt.plot(f_0,Pxx_den_0)
plt.plot(f_2,Pxx_den_2)
plt.plot(f_2no,Pxx_den_2no)
plt.plot(f_4,Pxx_den_4)
plt.plot(f_4no,Pxx_den_4no)
plt.grid()
plt.xlabel('cycle/month')
plt.ylabel('Power spectrum')
plt.legend(['periodogram','welch hanning default','welch segm=2, 50%overlap', \
           'welch segm=2, no overlap','welch segm=4, 50%overlap', \
           'welch segm=4, no overlap'])
plt.title('NW -- 20:30S - 105:115E')
plt.xlim([0, 0.2])
plt.show()


