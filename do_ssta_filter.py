'''
   calc the ssta on the filtered (spatially) data
'''

# Load required modules
import numpy as np
from scipy import signal
import pandas as pd
import xarray as xr

import datetime
from dateutil.relativedelta import relativedelta
from calendar import monthrange

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo

outfile_sstad = '/home/ecougnon/ana/SSTa_daily_filter_Aus_20032017.nc'
outfile_sstam = '/home/ecougnon/ana/SSTa_monthly_extremes_filter_Aus_20032017.nc'
fname = '/home/ecougnon/ana/SST_smooth_1deg_Aus.nc'

###########
# load AVHRR OISST smoothed
###########
ds_oisst = xr.open_dataset(fname)['__xarray_dataarray_variable__']. \
              sel(time=slice('2003-01-01','2017-10-31'))
#                  lat=slice(-50,-45), lon=slice(90,95))

# time period
time_vec = pd.date_range('2003-01-01','2017-10-31',name='time')
month_vec = pd.date_range('2003-01-01','2017-10-31',name='time',freq='M')

###########################################
# allocate memory fo rthe main variable
############################################
X = len(ds_oisst.lon)
Y = len(ds_oisst.lat)
nm = len(month_vec)
nd = len(time_vec)
'''
   Variable description from SST

    TMM -- Temperature monthly mean
    TMX -- Temperature monthly max
    TMN -- Temperature monthly min
    Tp90 -- Temperature monthly 90th percentile
    Tp10 -- Temperature monthly 10th percentile
'''
SSTa_d = xr.Dataset({'SSTa':(('time','lat','lon'),np.zeros((nd, Y, X)))}, \
                    {'time': time_vec, 'lat': ds_oisst.lat, 'lon': ds_oisst.lon})

SSTa_m = xr.Dataset({'TMM':(('time','lat','lon'),np.zeros((nm, Y, X))), \
                     'Tp90':(('time','lat','lon'),np.zeros((nm, Y, X))), \
                     'Tp10':(('time','lat','lon'),np.zeros((nm, Y, X)))}, \
                    {'time': month_vec, 'lat': ds_oisst.lat, 'lon': ds_oisst.lon})

################################################
# daily SSTa -- remove mean, trend and 
#               seasonal cycle
################################################
# deseason
dsea, sea, beta = np.apply_along_axis(eo.deseason_harmonic,0, \
                                      ds_oisst[:,:,:],4,365)
# detrend the time series
for ii in range(0,Y):
    for jj in range(0,X):
        valid = ~np.isnan(dsea[ii,jj])
        if (valid.any()==True):
            SSTa_d['SSTa'][:,ii,jj] = signal.detrend(np.squeeze(dsea[ii,jj]))
        elif (valid.all()==False):
            SSTa_d['SSTa'][:,ii,jj] = np.nan

SSTa_d.to_netcdf(outfile_sstad)

################################################
# monthly SSTa extremes/stats 
################################################
dtime = [datetime.date(2003,1,1) + datetime.timedelta(days=i) \
         for i in range(nd)]
str_id=np.empty(nm)
end_id=np.empty(nm)
k=0
for yy in range(2003,2017+1):
    if yy<2017:
        for mm in range(1,12+1):
            tmp1 = np.where(dtime == np.datetime64(datetime.date(yy,mm,1)))
            str_id[k] = tmp1[0][0]
            tmp2 = np.where(dtime == np.datetime64(datetime.\
                                                   date(yy,mm, \
                                                        monthrange(yy,mm)[1])))
            end_id[k] = tmp2[0][0]
            k = k +1
    elif yy==2017:
        for mm in range(1,10+1):
            tmp1 = np.where(dtime == np.datetime64(datetime.date(yy,mm,1)))
            str_id[k] = tmp1[0][0]
            tmp2 = np.where(dtime == np.datetime64(datetime.\
                                                   date(yy,mm, \
                                                        monthrange(yy,mm)[1])))
            end_id[k] = tmp2[0][0]
            k = k +1

SSTa_m['TMM'][:,:,:] = SSTa_d['SSTa'].resample('1MS', dim='time', how='mean')
for mm in range(0,nm):
    SSTa_m['Tp10'][mm,:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                                 SSTa_d['SSTa'][int(str_id[mm]): \
                                                                int(end_id[mm]+1), \
                                                                :,:],10)
    SSTa_m['Tp90'][mm,:,:] = np.apply_along_axis(np.nanpercentile,0, \
                                                 SSTa_d['SSTa'][int(str_id[mm]): \
                                                                int(end_id[mm]+1), \
                                                                :,:],90)

SSTa_m.to_netcdf(outfile_sstam)















