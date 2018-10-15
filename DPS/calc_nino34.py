'''
	calc the NINO index 3.4 (5degS-5degN ; 170-120degW or 190-240degE)
	from the DPS (enkf-9 refered to the climatology 
	based on 2003-2017) and the observation for the same based period

	Maybe do on enkf-13 -- forecast run for the assimilated
	period (2003-2006 ; 2002 is excluded) and the forecast 
	period (2007-2010) -- but not sure if very useful! For now
	focus on the reanalysis run (enkf-9)

	NINO3.4 is based on the SST anomalie from the above region.
	The anomalies are usually computed relative to a base period
	of 30years (problem in our DPS study! -- do with 2003-2017 on the enkf-9)

	Created: Apr 2018
	Author: Eva C.
	Last modified
'''


# Load required modules

import numpy as np
from scipy import io # read/load mat files
from scipy import signal  # detrend
from scipy.stats import pearsonr
import xarray as xr
import pandas as pd

from matplotlib import pyplot as plt
import datetime

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo
import eac_useful as eac

############################################
# reanalysis enkf-9
############################################

# load data
header_out = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'
#outfile = header_out + 'nino34_enkf9_mem001_20032017.nc'
# define the region 
lat_min = -5
lat_max = 5
# !!! LONGITUDE in the model!!! from -280 to 80 deg E
# NINO4: 160-210E
# NINO3.4: 190-240E
lon_min = 190 - 360 #160 -360 #190 - 360 #
lon_max = 240 - 360 #210 - 360 #240 - 360 #
yt_ocean = xr.open_dataset(gname)['yt_ocean']
yt_ocean = yt_ocean.sel(yt_ocean=slice(lat_min,lat_max))
xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
xt_ocean = xt_ocean.sel(xt_ocean=slice(lon_min,lon_max))

# monthly mean for the NINO34 region
SST_m = eac.read_netcdfs('/home/ecougnon/data/DPS/reanalysis/ETKF/ocean_daily_SST_enkf9_mem001_20??-20??.nc', \
                         dim='time', lat=yt_ocean, lon=xt_ocean, \
                         transform_func=lambda ds: ds.resample('1MS', \
                                                               dim='time', \
                                                               how='mean'))
SST_m = SST_m.squeeze('st_ocean')
SST_m = SST_m.sel(time=slice('2003-01-01','2017-11-30'))
SST_m = SST_m.mean(dim=('yt_ocean','xt_ocean'))

time_ts = pd.date_range('2003-01-01','2017-11-30',name='time',freq='M')

# calc monthly climatology
clim_m = SST_m.groupby('time.month').mean('time')
ssta_m = SST_m.groupby('time.month') - clim_m
#dssta_m = np.apply_along_axis(signal.detrend,0,ssta_m['temp'], type='linear')

'''
####################################################################
# compare with observations to evaluate the value of NINO3.4 against obs
# using HadISST as they're already monthly and with similar horizontal
# resolution (1deg)
####################################################################
ds_hadisst = xr.open_dataset('/home/data/sst/HadISST_sst.nc')['sst']. \
                sel(latitude=slice(lat_max,lat_min), \
                    longitude=slice(lon_min,lon_max), \
                    time=slice('2003-01-01','2017-11-30'))
ds_hadisst = ds_hadisst.mean(dim=('latitude','longitude'))
clim_obs = ds_hadisst.groupby('time.month').mean('time')
ssta_obs = ds_hadisst.groupby('time.month') - clim_obs
dssta_obs = np.apply_along_axis(signal.detrend,0,ssta_obs, type='linear')

#r1, _ = pearsonr(np.mean(dssta_m,axis=(1,2)),np.mean(dssta_obs,axis=(1,2)))

#linregress(np.mean(dssta_m,axis=(1,2)),np.mean(dssta_obs,axis=(1,2)))
'''
#'''
######################################
# compare with the smoothed (1/4 degree to 1degree)
# NOAA OISSTV2 daily dataset
#########################################
ds_oisst = xr.open_dataset('/home/ecougnon/ana/SST_smooth_1deg_NINO34region.nc')['__xarray_dataarray_variable__']. \
              sel(lat=slice(lat_min,lat_max), \
                  lon=slice(lon_min+360,lon_max+360), \
                  time=slice('2003-01-01','2017-11-30')). \
              resample('1MS', dim='time', how='mean')
ds_oisst = ds_oisst.mean(dim=('lat','lon'))
clim_oisst = ds_oisst.groupby('time.month').mean('time')
ssta_oisst = ds_oisst.groupby('time.month') - clim_oisst
#dssta_oisst = np.apply_along_axis(signal.detrend,0,ssta_oisst, type='linear')

#r2, _ = pearsonr(np.mean(dssta_m,axis=(1,2)),np.mean(dssta_oisst,axis=(1,2)))

#linregress(np.mean(dssta_m,axis=(1,2)),np.mean(dssta_oisst,axis=(1,2)))
#'''

############################################
# forecast enkf-13
############################################
fname_for = '/home/ecougnon/data/DPS/forecast/NINO34/sst_daily_PacEq_NINO_ensemble.nc'
sst_for = xr.open_dataset(fname_for)['sst'].sel(yt_ocean=slice(lat_min,lat_max), \
                                                xt_ocean=slice(lon_min,lon_max), \
                                                time=slice('2007-01-01', \
                                                           '2010-12-31')). \
                                            resample('1MS',dim='time',how='mean')
sst_for = sst_for.mean(dim=('yt_ocean','xt_ocean'))

ssta_for_enkf9 = sst_for.groupby('time.month') - clim_m
#dssta_for_enkf9 = np.apply_along_axis(signal.detrend,0,ssta_for_enkf9['temp'], \
#                                      type='linear')
'''
ssta_for_obs = sst_for.groupby('time.month') - clim_oisst
dssta_for_obs = np.apply_along_axis(signal.detrend,0,ssta_for_obs, \
                                    type='linear')

clim_for = sst_for.groupby('time.month').mean('time')
ssta_for = sst_for.groupby('time.month') - clim_for
dssta_for = np.apply_along_axis(signal.detrend,0,ssta_for, type='linear')
'''

###############################################
# plotting
##############################################
plt.figure(figsize=(11,4))
ax = plt.gca()
#plt.plot(time_ts,dssta_obs,'g')
plt.plot(time_ts,ssta_oisst,'k') #ds_oisst,'k') #
plt.plot(time_ts,ssta_m['temp'],'b') #SST_m['temp'],'b') #
plt.plot(sst_for.time,ssta_for_enkf9['temp'],'--b') #sst_for,'--b')#
#plt.plot(sst_for.time,dssta_for_obs,'.b')
#plt.plot(sst_for.time,sst_for,'b') #dssta_for,'-b')
plt.legend(['AVHRR OISST V2 (filtered to 1deg)',\
#            'NINO4 HadISST (1deg)', \
            'enkf-9', \
#            'NINO34 AVHRR OISST V2 (filtered to 1deg)',\
#            'NINO4 enkf13', \
            'forecast with clim. from enkf9'], loc=4) #\
#            'forecast with clim from Hadisst', \
#            'forecast with its own clim (2007-2010)'],loc=1)
#plt.legend(['NINO34 enkf-9',  \
#            'NINO34 AVHRR OISST V2 (filtered to 1deg)', \
#            'forecast with clim. from oisst (2003-2010)'])
plt.xlim(np.min(time_ts), np.max(time_ts))
#plt.ylim(-3,4)
#ax.annotate('r_hadisst = {:.2f}'.format(r1), xy=(.1,.09), xycoords=ax.transAxes)
#ax.annotate('r_avhrroisst = {:.2f}'.format(r2), xy=(.1,.03), xycoords=ax.transAxes)
'''
ax.fill_between(time_ts, ssta_m['temp'], \
                where=dssta_m >= 0.4, facecolor='red', \
                alpha=0.3, interpolate=True)
ax.fill_between(time_ts, ssta_m['temp'], \
                where=dssta_m <= -0.4, facecolor='blue', \
                alpha=0.3, interpolate=True)
'''
plt.title('SSTa from the NINO3.4 region')
plt.grid()
figfile = '/home/ecougnon/data/DPS/reanalysis/ETKF/nino34_SSTa.png'
#plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

'''
###################
# comparison graph
####################
plt.figure(figsize=(7,5))
plt.scatter(dssta_oisst[:-2],dssta_for, s=30, \
            c=np.arange(1,46+1), marker='*', \
            cmap=plt.cm.viridis)
#plt.scatter(dssta_m[48:94],dssta_for, s=30, c=np.arange(1,46+1), marker='o', \
#            cmap=plt.cm.viridis, facecolors=None)
plt.scatter(dssta_m[:-2],dssta_for, s=30, c=np.arange(1,46+1), marker='o', \
            cmap=plt.cm.viridis, facecolors=None)
cb=plt.colorbar()
cb.set_label('# of months')
plt.xlim([-3, 3])
plt.xlabel('enkf13 and oisst NINO3.4 (2007-2010 clim.)')
plt.ylim([-3, 3])
plt.ylabel('nino3.4 from the forecast run enkf13')
plt.legend(['enkf13 VS enkf9','enkf13 VS oisst'])
plt.title('NINO3.4 index')
plt.grid()

figfile_ = '/home/ecougnon/data/DPS/reanalysis/ETKF/nino34_2007_2010_clim_all_comp.png'
plt.savefig(figfile_, bbox_inches='tight', format='png', dpi=300)
'''

plt.show()






