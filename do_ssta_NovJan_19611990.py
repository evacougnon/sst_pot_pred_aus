'''
    calculate/test the climatology of the OISST V2 dataset by removing 
    the 1982-2005 mean with an offset calculated from the HadISST
    dataset by doing the difference between the SST mean from
    1961-1990 and 1982-2005

    USING only Nov-Jan months to be used in Sarah PK paper

    Author: Eva Cougnon
    Date created: Mar 2018
    Modified: 
'''
# load required modules

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
import eac_useful as eac


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

####################################################
## define the offset by using the HadISST dataset
## just for Nov to Jan of each years
####################################################
ds_hadisst_all = xr.open_dataset('/home/data/sst/HadISST_sst.nc')['sst']. \
                    sel(latitude=slice(lat_px_max+1,lat_px_min-1), \
                        longitude=slice(lon_px_min-1,lon_px_max+1), \
                        time=slice('1871-01-01','2018-05-31'))
lon = ds_hadisst_all.longitude
lat = ds_hadisst_all.latitude
# get the distance between each grid points in meters    
dx, dy = eac.dxdy(lat,lon)

#################################
# using latlon2km by Eric Oliver
#################################
'''
lon_test = ds_hadisst_all.longitude
lat_test = ds_hadisst_all.latitude
dx_test, dy_test = eo.dxdy(lon_test,lat_test) # in meters
'''
ds_hadisst = xr.open_dataset('/home/data/sst/HadISST_sst.nc')['sst']. \
                sel(latitude=slice(lat_px_max,lat_px_min), \
                    longitude=slice(lon_px_min,lon_px_max), \
                    time=slice('1871-01-01','2018-05-31'))#. \
#                mean(dim=('longitude','latitude'))
#'''
mask_land = np.ones((len(ds_hadisst.latitude),len(ds_hadisst.longitude)))
tmp = np.squeeze(ds_hadisst[0,:,:])
mask_land = np.ma.masked_where(np.isnan(tmp)==True, mask_land)

area = dx * dy * mask_land
#area_test = dx_test[:-1,:-1] * dy_test[:-1,:-1] * mask_land

# if working with Eric's dxdy use area_test instead of area
# area is based on WGS84 ellipsoid (geopy library
weight = area/np.nansum(area) 
tot_weight = weight.sum()
#'''
ds_hadisst = ((ds_hadisst*weight).sum(dim=('longitude','latitude'))) 
#ds_hadisst_clim = ds_hadisst.sel(time=slice('1961-01-01','1990-12-31'))

# average for each Nov-Jan of each year
# The first year will include only Jan and the last for this time series
# finishing at Dec 2017 will have only Nov-Dec
month_ = ds_hadisst.groupby('time.month').apply(lambda x: x).month
summer = (month_ == 11) | (month_ == 12) | (month_ == 1)
ds_hadisst_mean_NovJan = ds_hadisst.where(summer).resample('AS-Nov', 'time')

ds_hadisst_19611990 = ds_hadisst_mean_NovJan.sel(time=slice('1961-11-01', \
                                                            '1990-11-01')). \
                                             mean(dim=('time'))
ds_hadisst_19822005 = ds_hadisst_mean_NovJan.sel(time=slice('1982-11-01', \
                                                            '2005-11-01')). \
                                             mean(dim=('time'))
ds_hadisst_offset = ds_hadisst_19822005 - ds_hadisst_19611990

ds_hadisst_ssta = ds_hadisst_mean_NovJan - ds_hadisst_19611990

print('Offset between the 3 months average (Nov-Jan) of the period 1961-1990 and 1982-2005: ', ds_hadisst_offset, 'oC') 

#############################################
# extracting the SST data from AVHRR OISST V2
# run extract_region_sst_timeseries.py
#############################################
header = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/'
ds_oisst_d = xr.open_dataset(header + 'TasSea_Box_AreaAve_updated.nc')['sst']
ds_oisst_m = ds_oisst_d.resample('1MS', dim='time', how='mean')

month_m = ds_oisst_m.groupby('time.month').apply(lambda x: x).month
summer_m = (month_m == 11) | (month_m == 12) | (month_m == 1)
ds_oisst_m_annmean_NovJan = ds_oisst_m.where(summer_m).resample('AS-Nov', 'time')

month_d = ds_oisst_d.groupby('time.month').apply(lambda x: x).month
summer_d = (month_d == 11) | (month_d == 12) | (month_d == 1)
ds_oisst_d_annmean_NovJan = ds_oisst_d.where(summer_d).resample('AS-Nov', 'time')


ds_oisst_m_NovJan_19822005 = ds_oisst_m_annmean_NovJan.sel(time=slice('1982-11-01', \
                                                                      '2005-11-01')). \
                                                       mean(dim=('time'))

ds_oisst_d_NovJan_19822005 = ds_oisst_d_annmean_NovJan.sel(time=slice('1982-11-01', \
                                                                      '2005-11-01')). \
                                                       mean(dim=('time'))

ssta_d = ds_oisst_d_annmean_NovJan-(ds_oisst_d_NovJan_19822005-ds_hadisst_offset)
ssta_m = ds_oisst_m_annmean_NovJan-(ds_oisst_m_NovJan_19822005-ds_hadisst_offset)


print('SSTa from the daily dataset for Nov2017-Jan2018 (ds_oisst_d_annmean_NovJan[-1] - (ds_oisst_d_NovJan_19822005-ds_hadisst_offset)) ', ds_oisst_d_annmean_NovJan[-1] - (ds_oisst_d_NovJan_19822005-ds_hadisst_offset))

print('SSTa from the monthly dataset for Nov2017-Jan2018 (ds_oisst_m_annmean_NovJan[-1] - (ds_oisst_m_NovJan_19822005 - ds_hadisst_offset)) ', ds_oisst_m_annmean_NovJan[-1] - (ds_oisst_m_NovJan_19822005 - ds_hadisst_offset))

print('SSTa from the daily dataset for Nov2015-Jan2016 (ds_oisst_d_annmean_NovJan[-3] - (ds_oisst_d_NovJan_19822005-ds_hadisst_offset)) ', ds_oisst_d_annmean_NovJan[-3] - (ds_oisst_d_NovJan_19822005 - ds_hadisst_offset))

print('SSTa from the monthly dataset for Nov2015-Jan2016 (ds_oisst_m_annmean_NovJan[-3] - (ds_oisst_d_NovJan_19822005 - ds_hadisst_offset)) ', ds_oisst_m_annmean_NovJan[-3] - (ds_oisst_m_NovJan_19822005 - ds_hadisst_offset))

'''
######################################
# Saving in a netcdf file for Sarah PK to use
######################################
ssta_aa_novjan = xr.Dataset({'SST_d':(['time_d'], ds_oisst_d_annmean_NovJan), \
                             'SSTa_d':(['time_d'], ssta_d), \
                             'SST_m':(['time_m'], ds_oisst_m_annmean_NovJan), \
                             'SSTa_m':(['time_m'], ssta_m), \
                             'offset':(ds_hadisst_offset)}, \
                            coords={'time_d': ds_oisst_d_annmean_NovJan.time, \
                                    'time_m': ds_oisst_m_annmean_NovJan.time})
ssta_aa_novjan.attrs['title'] = 'Yearly NovJan averaged (1982-2018) icalculated with the daily and monthly OISSTV2, including the offset to calculate the climatology based on 1961-1990, for the Tasman Sea Box (135oE). Area averaged using WGS84 to calculate the distances between the regular 1/4 degree grid'

#ssta_sa_novjan.attrs['title'] = 'Yearly NovJan averaged (1982-2018) SST and SSTa with the daily and monthly OISSTV2, including the offset to calculate the climatology based on 1961-1990, for the Tasman Sea Box (135oE). Spatially averaged considering a 1/4 degree regular grid'

ssta_aa_novjan.to_netcdf(header + 'ssta_aa_novjan_TasSeaBox135_updated.nc')

'''
############################################
# plotting the NovJan anomalies
############################################
plt.figure(figsize=(13,3))
ax = plt.subplot(111)
plt.plot(np.arange(1982,2017+1),ssta_d[1:],'-*b')
plt.plot(np.arange(1871,2018), ds_hadisst_ssta[1:],'-*k')
#plt.plot(ds_oisst_d_annmean_NovJan.time[1:],ssta_d[1:],'*k')
#ax.set_xlim(['1982-01-01','2018-03-25'])
ax.set_xlim([1871,2018])
ax.set_ylim([-2, 2])
plt.legend(['NOAA OISST v2', 'HadISST'])
plt.title('NovJan SSTa based on 1961-1990 climatology')
#plt.title('NovJan SSTa based on the 1982-2005 climatology adjusted to be relative to 1961-1990', fontsize=16)
plt.grid()

#figfile_eps = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SSTa_NovJan_19822017_updated.eps'
#plt.savefig(figfile_eps, bbox_inches='tight', format='eps', dpi=300)

figfile = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SSTa_NovJan_NOAAOISST_HadISST.png'
plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

plt.show()



