'''

    Adapted from Eric's code gen_annual_extremes
    for monthly extreme calculated from daily data

'''

# Load required modules

import numpy as np
from scipy import io # read/load mat files
from scipy import linalg
from scipy import stats
import datetime
import time as time

from dateutil.relativedelta import relativedelta
from calendar import monthrange
from netCDF4 import Dataset

import eric_oliver as eo

# Load some select time series
#
# observations
#
pathroot = '/home/ecoliver/'
#pathroot = '/home/ecoliver/Desktop/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
file0 = header + '1982/avhrr-only-v2.19820101.nc'
t = np.arange(datetime.date(1982,1,1).toordinal(), \

T = len(t)
year = np.zeros((T))
month = np.zeros((T))
day = np.zeros((T))
for i in range(T):
    year[i] = datetime.date.fromordinal(t[i]).year
    month[i] = datetime.date.fromordinal(t[i]).month
    day[i] = datetime.date.fromordinal(t[i]).day
###
# time period
MinYear = 1982
MaxYear = 2016
NumYears = MaxYear - MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
dtime = [datetime.date(MinYear,1,1) + datetime.timedelta(days=i) \
         for i in range(NumDays)]
years = np.unique(year)
months = np.unique(month)
NB = len(years)*len(months)

# start and end index for each month
str_id=np.empty(NB)
end_id=np.empty(NB)
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
str_id=str_id.astype(int)
end_id=end_id.astype(int)
#
# lat and lons of obs
#
fileobj = Dataset(file0, mode='r')
lon_tmp = fileobj.variables['lon'][:].astype(float)
lat_tmp = fileobj.variables['lat'][:].astype(float)
# not use
#fill_value = fileobj.variables['sst']._FillValue.astype(float)
#scale = fileobj.variables['sst'].scale_factor.astype(float)
#offset = fileobj.variables['sst'].add_offset.astype(float)
fileobj.close()

#
# useful numbers
res = 0.25 # resolution of the data -- quarter degree
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
# allocate memory fo rthe main variable
X = len(lon_map)
Y = len(lat_map)
SST = {}
keys = ['TMM', 'TMX', 'TMN', 'Tp90', 'Tp10'] 
'''
   Variable description from SST

    TMM -- Temperature monthly mean
    TMX -- Temperature monthly max
    TMN -- Temperature monthly min
    Tp90 -- Temperature monthly 90th percentile
    Tp10 -- Temperature monthly 10th percentile
'''
for key in keys:
    SST[key] = np.zeros((Y, X, NB))

#
# loop through locations
#
icnt = 0
for i in lon_id:
    t_lon = time.time()
    print(icnt+1, 'of', len(lon_map))
#   load SST
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + \
                        str(i+1).zfill(4) + '.mat')
    sst_ts = matobj['sst_ts']
!!!    sst_ts = np.array(sst_ts) !! may want to change the dimension 
!!!    timexlat instead of lat by time
#    sst_ma = np.ma.masked_where(np.isnan(sst_ts) == True, sst_ts, copy=True)
#    jcnt = 0
#    for j in lat_id:
    tmp = np.apply_along_axis(eo.deseason_harmonic,1,sst_ts[lat_id,:],4,365)
#tmp is a matrix (lat,3) 3 are the 3 output from the function with a length of time
  dsea_sst = np.array(tmp[:,0][:])
# incomplete but that's the way to go!!
    for mm in range(0,NB):
         SST['TMM'][jcnt,icnt,mm] = np.mean(dsea_sst[int(str_id[mm]) \
                                                     :int(end_id[mm]+1)])
         SST['TMX'][jcnt,icnt,mm] = np.max(dsea_sst[int(str_id[mm]) \
                                                    :int(end_id[mm]+1)])
         SST['TMN'][jcnt,icnt,mm] = np.min(dsea_sst[int(str_id[mm]) \
                                                    :int(end_id[mm]+1)])
         SST['Tp90'][jcnt,icnt,mm] = np.percentile(dsea_sst[int(str_id[mm]) \
                                                            :int(end_id[mm]+1)] \
                                                   , 90)
         SST['Tp10'][jcnt,icnt,mm] = np.percentile(dsea_sst[int(str_id[mm]) \
                                                            :int(end_id[mm]+1)] \
                                                   , 10)
        jcnt += 1
    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt += 1

'''
    jcnt = 0
    for j in lat_id:
# checking that the time series is not empty -- deseason does not
# handle a complete nan vector
        if ((np.isnan(np.min(sst_ts[j,:])) == False) & \
            (np.isnan(np.max(sst_ts[j,:])) == False)):

            dsea_sst,season,beta = eo.deseason_harmonic(sst_ts[j,:],4,365)
        
            for mm in range(0,NB):


                SST['TMM'][jcnt,icnt,mm] = np.mean(dsea_sst[int(str_id[mm]) \
                                                            :int(end_id[mm]+1)])
                SST['TMX'][jcnt,icnt,mm] = np.max(dsea_sst[int(str_id[mm]) \
                                                           :int(end_id[mm]+1)])
                SST['TMN'][jcnt,icnt,mm] = np.min(dsea_sst[int(str_id[mm]) \
                                                           :int(end_id[mm]+1)])
                SST['Tp90'][jcnt,icnt,mm] = np.percentile(dsea_sst[int(str_id[mm]) \
                                                                   :int(end_id[mm]+1)] \
                                                          , 90)
                SST['Tp10'][jcnt,icnt,mm] = np.percentile(dsea_sst[int(str_id[mm]) \
                                                                   :int(end_id[mm]+1)] \
                                                          , 10)
        jcnt += 1
    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt += 1
'''

# Save data 
#outfile = '/home/ecougnon/ana/SSTa_monthly_extremes_Aus'
#np.savez(outfile, lon_map=lon_map, lat_map=lat_map, SST=SST)

####
# NOTES
# when looping through lon/lat and time after removing the seasons, 
# each loop through lon takes > 70s 
#
# when using a masked array > 75s
#



