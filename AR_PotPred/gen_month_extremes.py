'''

    Adapted from Eric's code gen_annual_extremes
    for monthly extreme calculated from daily data

'''

# Load required modules

import numpy as np
from scipy import io # read/load mat files
import pandas as pd
import xarray as xr
import time as time
import datetime
from dateutil.relativedelta import relativedelta
from calendar import monthrange

from netCDF4 import Dataset

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

t_ini = time.time()
# saving file -- with nc, does not allow to overwrite!!!
outfile = '/home/ecougnon/ana/SST_monthly_extremes_Aus.nc'
# time period
time_vec = pd.date_range('1982-01-01','2016-12-31',name='time')
month_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
NumYears = 2016-1982+1
# time index to calc percentile
MinYear=1982
MaxYear=2016
NumYears = MaxYear - MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
dtime = [datetime.date(MinYear,1,1) + datetime.timedelta(days=i) \
         for i in range(NumDays)]
str_id=np.empty(len(month_vec))
end_id=np.empty(len(month_vec))
k=0
for yy in range(MinYear,MaxYear+1):
    for mm in range(1,12+1):
        tmp1 = np.where(dtime == np.datetime64(datetime.date(yy,mm,1)))
        str_id[k] = tmp1[0][0]
        tmp2 = np.where(dtime == np.datetime64(datetime.\
                                               date(yy,mm, \
                                                    monthrange(yy,mm)[1])))
        end_id[k] = tmp2[0][0]
        k = k +1
#
# choose lats and lons
#
pathroot = '/home/ecoliver/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
file0 = header + '1982/avhrr-only-v2.19820101.nc'
fileobj = Dataset(file0, mode='r')
lon_tmp = fileobj.variables['lon'][:].astype(float)
lat_tmp = fileobj.variables['lat'][:].astype(float)
fileobj.close()
# resolution of the data -- quarter degree
res = 0.25 
# define region of study 
lat_min = -55 #-17 #-55
lat_max = 10.25 #0 #10.25
lon_min = 90 #155 #90
lon_max = 180.25 #165 #180.25
# find the index of the lat/lon min/max
lat_min_id=np.nonzero((lat_tmp>(lat_min-res)) & (lat_tmp<(lat_min)))
lat_max_id=np.nonzero((lat_tmp>(lat_max-res)) & (lat_tmp<(lat_max)))
lon_min_id=np.nonzero((lon_tmp>(lon_min-res)) & (lon_tmp<(lon_min)))
lon_max_id=np.nonzero((lon_tmp>(lon_max-res)) & (lon_tmp<(lon_max)))
# get the matrix of indexes for the lat/lon
lat_map = lat_tmp[lat_min_id[0][0]:lat_max_id[0][0]+1]
lat_id = range(lat_min_id[0][0],lat_max_id[0][0]+1,1)
lon_map = lon_tmp[lon_min_id[0][0]:lon_max_id[0][0]+1]
lon_id = range(lon_min_id[0][0],lon_max_id[0][0]+1,1)

#
# allocate memory fo rthe main variable
#
X = len(lon_map)
Y = len(lat_map)
nm = len(month_vec) 
#SST = {}
#keys = ['TMM', 'TMX', 'TMN', 'Tp90', 'Tp10'] 
'''
   Variable description from SST

    TMM -- Temperature monthly mean
    TMX -- Temperature monthly max
    TMN -- Temperature monthly min
    Tp90 -- Temperature monthly 90th percentile
    Tp10 -- Temperature monthly 10th percentile
'''
#for key in keys:
#    SST[key] = np.zeros((Y, X, t))
SST_ = xr.Dataset({'TMM':(('time','lat','lon'),np.zeros((nm, Y, X))), \
                   'TMX':(('time','lat','lon'),np.zeros((nm, Y, X))), \
                   'TMN':(('time','lat','lon'),np.zeros((nm, Y, X))), \
                   'Tp90':(('time','lat','lon'),np.zeros((nm, Y, X))), \
                   'Tp10':(('time','lat','lon'),np.zeros((nm, Y, X)))}, \
                  {'time': month_vec, 'lat': lat_map, 'lon': lon_map})
elapsed_ini = time.time() - t_ini
print('elapsed time for initialisation:', elapsed_ini)
#
# loop through longitude
#
icnt = 0
for i in lon_id:
    t_lon = time.time()
    print(icnt+1, 'of', len(lon_map))
#   load SST
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + \
                        str(i+1).zfill(4) + '.mat')
    
    ds_sst = xr.Dataset({'sst_ts':(('lat','time'), \
                                   matobj['sst_ts'][lat_id,:])}, \
                        {'time': time_vec,'lat': matobj['lat'][lat_id,0]})
    
    monthly_avg_sst = ds_sst.resample('1MS', dim='time', how='mean')
    monthly_max_sst = ds_sst.resample('1MS', dim='time', how='max')
    monthly_min_sst = ds_sst.resample('1MS', dim='time', how='min')
#    monthly_p10_sst = ds_sst.resample('1MS', dim='time', \
#                                      how= lambda x: x.quantile(0.1))
####
# NOTES -- when doing only mean/max/min -- <6s per longitude
#          when including the following part with p10 and p90 -- >50s per lon
###
    SST_['TMM'][:,:,icnt]=monthly_avg_sst['sst_ts']
    SST_['TMX'][:,:,icnt]=monthly_max_sst['sst_ts']
    SST_['TMN'][:,:,icnt]=monthly_min_sst['sst_ts']

    for mm in range(0,len(month_vec)):
        SST_['Tp10'][mm,:,icnt] = ds_sst['sst_ts'][:,int(str_id[mm]): \
                                                     int(end_id[mm]+1)] \
                                  .quantile(0.1,dim='time')
        SST_['Tp90'][mm,:,icnt] = ds_sst['sst_ts'][:,int(str_id[mm]): \
                                                     int(end_id[mm]+1)] \
                                  .quantile(0.9,dim='time')

    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)

    icnt += 1

# Saving into a netcdf file
SST_.to_netcdf(outfile)


