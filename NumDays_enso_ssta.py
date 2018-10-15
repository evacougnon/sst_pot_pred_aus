# load SSTs daily time series and estimate a percentage of the data
# above/below p90/p10 values during El Nino and La Nina yeats



# load required modules

import numpy as np
from scipy import stats
import pandas as pd
import xarray as xr

from datetime import date
from calendar import monthrange
import time as time

# saving file
#outfile = '/home/ecougnon/ana/NumDays_ElN_p90p10_Aus.nc'
outfile = '/home/ecougnon/ana/NumDays_LaN_p90p10_Aus.nc'
# Load data
# Daily SSTa
fname = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
tim = xr.open_dataset(fname)['time']
# number of days for 10% of the total number of days
tim_10th = len(tim)/10
# SSTa time series
sst_ts = xr.open_dataset(fname)['SSTa']
# threshold values (from the daily SSTa time series)
fname2 = '/home/ecougnon/ana/SSTa_map_Aus.nc'
#SST_X = xr.open_dataset(fname2)
sst_p10 = np.array(xr.open_dataset(fname2)['Tp10'])
sst_p90 = np.array(xr.open_dataset(fname2)['Tp90'])

# usefull numbers
MinYear = 1982
MaxYear =2016
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
# time vector for plotting
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,31).toordinal()+1)
# allocate memory
Y = len(lon)
X = len(lat)
'''
n_p90_LaN -- number of days above p90 corresponding to La Nina year
n_p10_LaN -- number of days below p10 corresponding to La Nina year
p_p90_LaN -- percentage of number of days above p90 corresponding to La Nina year
p_p10_LaN -- percentage of number of days below p10 corresponding to La Nina year
'''
Pdays = xr.Dataset({'n_p90_LaN':(('lon','lat'),np.zeros((Y, X))), \
                    'n_p10_LaN':(('lon','lat'),np.zeros((Y, X))), \
                    'p_p90_LaN':(('lon','lat'),np.zeros((Y, X))), \
                    'p_p10_LaN':(('lon','lat'),np.zeros((Y, X)))}, \
                   {'lat': lat, 'lon': lon})
'''
Pdays = xr.Dataset({'n_p90_ElN':(('lon','lat'),np.zeros((Y, X))), \
                    'n_p10_ElN':(('lon','lat'),np.zeros((Y, X))), \
                    'p_p90_ElN':(('lon','lat'),np.zeros((Y, X))), \
                    'p_p10_ElN':(('lon','lat'),np.zeros((Y, X)))}, \
                   {'lat': lat, 'lon': lon})
'''
## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 

file_enso = np.genfromtxt('/home/ecougnon/data/enso/MEI_index.txt', \
                        skip_header=10, skip_footer = 30, delimiter='\t')
str_id = np.nonzero((file_enso[:,0]>(MinYear-1)) \
                    & (file_enso[:,0]<(MinYear+1)))

enso_id = np.zeros(NumYears*12) # the idexes to use for our study (from 1982)
dtime_enso = [None] * (NumYears*12)
enso_daily = np.empty(NumDays)
k=0
l=0
d=0
for yy in np.arange(str_id[0][0],len(file_enso[:,0])):
    for mm in np.arange(1,12+1):
        enso_id[k] = file_enso[yy,mm]
        dtime_enso[k] = date(MinYear+l,mm,1).toordinal()
        enso_daily[d:d+monthrange(MinYear+l,mm)[1]] = enso_id[k]
        k = k + 1
        d = d + monthrange(MinYear+l,mm)[1]
    l = l + 1
enso_p_id = np.nonzero(enso_daily>=0.75)
enso_n_id = np.nonzero(enso_daily<=-0.75)
enso_neu_id = np.nonzero((enso_daily>-0.75) & (enso_daily<0.75))

'''
# SSTa for positive ENSO
sst_ElN = np.array(sst_ts)
sst_ElN[:,:,enso_n_id] = np.nan
sst_ElN[:,:,enso_neu_id] = np.nan
'''
# SSTa for negative ENSO
sst_LaN = np.array(sst_ts)
sst_LaN[:,:,enso_p_id] = np.nan
sst_LaN[:,:,enso_neu_id] = np.nan


for ln in range(0,Y):
    for lt in range(0,X):
        
        Pdays['n_p90_LaN'][ln,lt] = np.count_nonzero(sst_LaN[ln,lt,:] \
                                                     > sst_p90[ln,lt], \
                                                     axis=2)
        Pdays['n_p10_LaN'][ln,lt] = np.count_nonzero(sst_LaN[ln,lt,:] \
                                                     < sst_p10[ln,lt], \
                                                     axis=2)
        Pdays['p_p90_LaN'][ln,lt] = Pdays['n_p90_LaN'][ln,lt]/tim_10th
        Pdays['p_p10_LaN'][ln,lt] = Pdays['n_p10_LaN'][ln,lt]/tim_10th 
        
        '''
        Pdays['n_p90_ElN'][ln,lt] = np.count_nonzero(sst_ElN[ln,lt,:] \
                                                     > sst_p90[ln,lt], \
                                                     axis=2)
        Pdays['n_p10_ElN'][ln,lt] = np.count_nonzero(sst_ElN[ln,lt,:] \
                                                     < sst_p10[ln,lt], \
                                                     axis=2)
        Pdays['p_p90_ElN'][ln,lt] = Pdays['n_p90_ElN'][ln,lt]/tim_10th
        Pdays['p_p10_ElN'][ln,lt] = Pdays['n_p10_ElN'][ln,lt]/tim_10th
        '''

## save files
Pdays.to_netcdf(outfile)



