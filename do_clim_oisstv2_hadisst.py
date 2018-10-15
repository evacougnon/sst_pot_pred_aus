'''
    calculate/test the climatology of the OISST V2 dataset by removing 
    the 1982-2005 mean with an offset calculated from the HadISST
    dataset by doing the difference between the SST mean from
    1961-1990 and 1982-2005

    this is following Eric's papers: Oiver et al 2017 Nat Com
    and Oliver et al. 2018 BAMS

    Author: Eva Cougnon
    Date created: Feb 2018
    Modified: Mar 2018

'''
# load required modules

import numpy as np
import xarray as xr
import pandas as pd
#import time as time
from datetime import date
from geopy.distance import vincenty # distance lat/lon

from scipy import io # load matlab file
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays
from scipy.stats import pearsonr

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm
import cmocean

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo
import marineHeatWaves as mhw
import eac_useful as eac

#figfile ='/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SSTa_oisst_hadisst_offset.png'

#####################
# define some conditions
##########################
# What to plot between ssta and mhw timeseries
plot = 'mhw' #'mhw' # ssta
# whether of not turing on the area-averaged (area-weighted)
# or just a spatial average 9no weight
what = 'S' #'A' # 'S'

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
lat_px_min = -45 #-45 #-46 # deg N
lat_px_max = -37 #-37 #-26
lon_px_min = 147 #139 #147 #135 #150 # deg E
lon_px_max = 155 #147 #155 #174

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

if plot=='ssta':
# get anomalie by removing a climatology based on 1961-1990
    dsea_tmp1, sea_tmp1, beta = eo.deseason_harmonic(np.array(ds_hadisst_clim), \
                                                     2,12)
    dsea_hadisst = np.empty(len(ds_hadisst_mean))
    dsea_hadisst.fill(np.nan)
    for tt in range(0,len(ds_hadisst_mean),12):
        for mm in range(0,12):
            dsea_hadisst[tt+mm] = np.array(ds_hadisst_mean[tt+mm]) - sea_tmp1[mm]
#############################################
# runnind standard deviation
#############################################
hadisst_std_30 = np.empty(len(ds_hadisst_mean))
hadisst_std_30.fill(np.nan)
#hadissta_std_30 = np.empty(len(dsea_hadisst))
#hadissta_std_30.fill(np.nan)
w = 30
for tt in range(0,len(tim_vec_hadisst)-int(w/2)):
    hadisst_std_30[tt+int(w/2)] = np.nanstd(ds_hadisst_mean[tt:tt+w])
#    hadissta_std_30[tt+int(w/2)] = np.nanstd(dsea_hadisst[tt:tt+w])
#############################################
# extracting the SST data from AVHRR OISST V2
# run extract_region_sst_timeseries.py
#############################################
header = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/'
if what == 'S':
    ds_oisst_d = xr.open_dataset(header + 'TasSea_135.nc')['sst']
elif what == 'A':
    ds_oisst_d = xr.open_dataset(header + \
                                 'TasSea_Box_AreaAve_updated_May.nc')['sst']

ds_oisst_m = ds_oisst_d.resample('1MS', dim='time', how='mean')
'''
 if using another region, make sure to use extract_region_sst_timeseries.py
'''

ds_oisst_19822005 = ds_oisst_m.sel(time=slice('1982-01-01','2005-12-31'))
tim_vec_oisst = pd.date_range('1982-01-01','2018-03-31',name='time',freq='M')
if plot=='ssta':
# get anomalie by removing a climatology based on 1982-2005
    dsea_tmp2, sea_tmp2, beta = eo.deseason_harmonic(np.array(ds_oisst_19822005), \
                                                     2,12)
    dsea_oisst = np.empty(len(ds_oisst_m))
    dsea_oisst.fill(np.nan)
    for tt in range(0,len(ds_oisst_m)-3,12):
        for mm in range(0,12):
            dsea_oisst[tt+mm] = np.array(ds_oisst_m[tt+mm]) - \
                                (sea_tmp2[mm] - np.array(ds_hadisst_offset))
# add Jan 2018
    dsea_oisst[-3] = np.array(ds_oisst_m[-3])-(sea_tmp2[0]-np.array(ds_hadisst_offset))
# add Feb 2018
    dsea_oisst[-2] = np.array(ds_oisst_m[-2])-(sea_tmp2[1]-np.array(ds_hadisst_offset))
# add Mar 2018
    dsea_oisst[-1] = np.array(ds_oisst_m[-1])-(sea_tmp2[2]-np.array(ds_hadisst_offset))

# without the offset 
    dsea_oisst_nooffset = np.empty(len(ds_oisst_m))
    dsea_oisst_nooffset.fill(np.nan)
    for tt in range(0,len(ds_oisst_m)-3,12):
        for mm in range(0,12):
            dsea_oisst_nooffset[tt+mm] = np.array(ds_oisst_m[tt+mm])-sea_tmp2[mm]

# add Jan 2018
    dsea_oisst_nooffset[-3] = np.array(ds_oisst_m[-3])-sea_tmp2[0]
# add Feb 2018
    dsea_oisst_nooffset[-2] = np.array(ds_oisst_m[-2])-sea_tmp2[1]
# add Mar 2018
    dsea_oisst_nooffset[-1] = np.array(ds_oisst_m[-1])-sea_tmp2[2]

elif plot=='mhw':
####################################
## apply the MHW code
####################################
#time vector for the mhw function!!!
# needs to be generated wth datetime format!
    MinYear = 1982
    MaxYear =2018 # warninng, so far finishes in Jan!!
    NumYears = MaxYear-MinYear+1
    MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                     # case the first year is a leap year
    NumDays = 365*NumYears + MaxNumLeapYear
# warning finishes in Feb 2018 not Dec
    dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,5,19).toordinal()+1)

    sst_mean_d = np.array(ds_oisst_d)
    dT = np.array(ds_hadisst_offset)
    mhws, clim = mhw.detect(dtime, sst_mean_d, climatologyPeriod=[1982,2005], \
                            alternateClimatology=[dtime, sst_mean_d - dT]) 

# time vector for plotting with mhw
    tim_vec_plot = pd.date_range('1982-01-01','2018-05-19',name='time',freq='D')

# SSTa mean from the daily data from 
    '''
    Ind = -146
    Ind_end = -29
    SSTa_area_mean = np.mean(sst_mean_d[Ind:Ind_end+1] - clim['seas'][Ind:Ind_end+1])
    print('SSTa area mean for the period starting from ' , tim_vec_plot[Ind], 'until' , tim_vec_plot[Ind_end], ': ')
    print('SSTa = ', SSTa_area_mean)
    '''
'''
extra info for indexes when using the file finishing in Feb 2018: 
1Nov2017: -120
15Nov2017: -106
1Dec2017: -90
1Jan2018: -59
31Jan2018: -29
1Feb2018: -28

1Oct2015: 12326
31Jul2016: 12630
'''


#####################################
## plotting
#####################################
if plot=='ssta':
    plt.figure(figsize=(13,18))
    ax = plt.subplot(311)
    plt.plot(ds_hadisst_mean.time,ds_hadisst_mean, 'k')
    plt.plot(ds_oisst_d.time,ds_oisst_d,'b')
    plt.legend(['monthly HadISST','daily OISST V2'])
    ax.set_xlim(['1871-01-01','2018-01-31'])
    ax.set_ylim([13, 21])
    plt.title('SST from HadISST and OISST V2')

    ax = plt.subplot(312)
#    plt.plot(tim_vec_plot, sst_mean_d-clim['seas'])
    plt.plot(tim_vec_oisst,dsea_oisst,'b')
    plt.plot(tim_vec_hadisst,dsea_hadisst, 'k')
    plt.legend(['monthly HadISST','daily OISST V2'])
    ax.set_xlim(['1871-01-01','2018-01-31'])
    ax.set_ylim([-1.5, 2.5])
    plt.title('SST anomalies (climatology based on 1961-1990)')

    ax=plt.subplot(313)
    plt.plot(tim_vec_hadisst,dsea_hadisst, 'k')
    plt.plot(tim_vec_oisst,dsea_oisst, 'b')
    plt.plot(tim_vec_oisst,dsea_oisst_nooffset, '--')
    plt.legend(['monthly HadISST','daily OISST V2 adjusted', \
                'daily OISST V2 no offset'])
    ax.set_xlim(['1982-01-01','2018-01-31'])
    ax.set_ylim([-1.5, 2.5])
    plt.title('SSTa', fontsize=16)
    tmp1=xr.Dataset({'dsea_oisst':(('time'),dsea_oisst)},{'time':tim_vec_oisst})
    tmp1=tmp1.resample('1MS',dim='time',how='mean')
    r1, _ = pearsonr(dsea_hadisst[1332:],tmp1['dsea_oisst'][:-3])
    tmp2=xr.Dataset({'dsea_oisst__nooffset':(('time'),dsea_oisst_nooffset)}, \
                    {'time':tim_vec_oisst})
    tmp2=tmp2.resample('1MS',dim='time',how='mean')
    r2, _ = pearsonr(dsea_hadisst[1332:],tmp2['dsea_oisst__nooffset'][:-3])
    ax.annotate('coef. corr. HadISST and NOAA OISST corrected = {:.2f}'.format(r1), \
                xy=(.1,.09), xycoords=ax.transAxes)
    ax.annotate('coef. corr. HadISST and NOAA OISST = {:.2f}'.format(r2), \
                xy=(.1,.03), xycoords=ax.transAxes)
#    plt.grid()
 
    plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)


    plt.figure(figsize=(13,4))
    ax = plt.subplot()
    plt.plot(tim_vec_hadisst, hadisst_std_30,'k')
    plt.title('SST 30-year running standard deviation')
    ax.set_xlim(['1900-01-01','2020-12-31'])
#    ax.set_ylim([1.7, 1.85])
    plt.grid()
    plt.show()

elif plot=='mhw':
# Find indices the shaded areas

    plt.figure(figsize=(13,18)) #13,3)) #13,18))
#    '''
    ax = plt.subplot(311)
    plt.plot(tim_vec_plot, clim['seas'], '0.5')
    plt.plot(tim_vec_plot, sst_mean_d, 'k')
    plt.plot(tim_vec_plot, clim['thresh'],'b')
    plt.legend(['climatology','SST daily','threshold'])
    ax.set_xlim(['2015-01-01','2018-05-31'])
    ax.set_ylim([13, 22])
    ax.fill_between(tim_vec_plot, sst_mean_d, \
                    where=sst_mean_d>= clim['thresh'], \
                     facecolor='red', alpha=0.5, \
                   interpolate=True)
    plt.title('SST Tasman Sea Oliver2017 area', fontsize=16)
    plt.grid()
    
    ax = plt.subplot(312)
    plt.plot(tim_vec_plot, sst_mean_d-clim['seas'], 'k')
    ax.set_xlim(['2015-01-01','2018-05-31'])
    ax.set_ylim([-1.5, 3.5])
    plt.axvline(x='2017-10-01')
    plt.title('SSTa Tasman Sea Oliver2017 area', fontsize=16)
    plt.grid()
#   
    ax = plt.subplot(313)
    plt.plot(tim_vec_plot, sst_mean_d-clim['seas'], 'k')
    ax.set_xlim(['2017-10-01','2018-05-31'])
    ax.set_ylim([-1.5, 3.5])
    plt.title('SSTa Tasman Sea Oliver2017, from blue line above till Feb2018', \
              fontsize=16)
    plt.grid()
    '''
    figfile = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SST_TasSeaBox_AreaAve_updated.png'
    plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

    figfile_eps = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/SST_TasSeaBox_AreaAve_updated.eps'
    plt.savefig(figfile, bbox_inches='tight', format='eps', dpi=300)
    '''
    plt.show()


############################################
# saving the data set for Sarah K.
# SST, SSTa, climatology time series
# averaged for the season
############################################
'''
tm_mhw = xr.Dataset({'SST':(['time'], sst_mean_d), \
                     'clim':(['time'], clim['seas']), \
                     'SSTa':(['time'], sst_mean_d-clim['seas'])}, \
                    coords={'time': tim_vec_plot})
tm_mhw.attrs['title'] = 'Timeseries (Jan1982-Feb2018) of SSTs, SSTa, climatology based on 1961-1990, for the Tasman Sea Box (135oE)'


tm_mhw.to_netcdf(header + 'tm_mhw_TasSeaBox135.nc')

'''

