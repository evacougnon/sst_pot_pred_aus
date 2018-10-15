'''
    compare weighted mean and non weighted mean of SST
'''



#load required modules

import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
from geopy.distance import vincenty # distance lat/lon

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


# define indexes for lat lon of one pixel
'''
info on the regions: 
        Tasman Sea Box with extended western boundary:
                lat[-46:-26]; lon[135:174]
        SE Aus, Oliver et al. 2017 box: 
                lat[-45:-37]; lon[147:155]
        west Tasmania, west of the above box:
                lat[-45:-37]; lon[139:147]
'''
lat_px_min = -46 #-45 #-46 # deg N
lat_px_max = -26 #-37 #-26
lon_px_min = 135 #139 #147 #135 #150 # deg E
lon_px_max = 174 #147 #155 #174


header = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/'

ds_sa = xr.open_dataset(header + 'TasSea_135.nc')['sst']
ds_aa = xr.open_dataset(header + 'TasSea_Box_AreaAve_test.nc')['sst']. \
           sel(time=slice('1982-01-01','2018-02-28'))

ds_ssta_sa = xr.open_dataset(header + 'ssta_sa_novjan_TasSeaBox135.nc')['SSTa_d']
ds_ssta_aa = xr.open_dataset(header + 'ssta_aa_novjan_TasSeaBox135.nc')['SSTa_d']

plt.figure()
plt.subplot(1,2,1)
plt.plot(ds_sa,ds_aa,'*')
plt.grid()
plt.xlabel('not-weighted SST')
plt.ylabel('weighted SST')
plt.xlim(13,22)
plt.ylim(10,19)
plt.title('daily SST for the region')

plt.subplot(1,2,2)
plt.plot(ds_ssta_sa[1:],ds_ssta_aa[1:],'*')
plt.grid()
plt.xlabel('not-weighted SSTa')
plt.ylabel('weighted SSTa')
plt.xlim(-0.5,1.8)
plt.ylim(-0.5,1.8)
plt.title('yearly Nov-Jan SSTa for the region')

plt.show()
 












