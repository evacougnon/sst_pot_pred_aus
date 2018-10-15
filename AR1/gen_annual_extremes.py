'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to 
  select SST time series around the globe

'''

# Load required modules

import numpy as np
from scipy import io
from scipy import linalg
from scipy import stats
from datetime import date
from netCDF4 import Dataset

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

# Load some select time series

#
# observations
#

pathroot = '/mnt/erebor/'
#pathroot = '/home/ecoliver/Desktop/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
file0 = header + '1982/avhrr-only-v2.19820101.nc'
t = np.arange(date(1982,1,1).toordinal(),date(2016,12,31).toordinal()+1)
T = len(t)
year = np.zeros((T))
month = np.zeros((T))
day = np.zeros((T))
for i in range(T):
    year[i] = date.fromordinal(t[i]).year
    month[i] = date.fromordinal(t[i]).month
    day[i] = date.fromordinal(t[i]).day

#
# lat and lons of obs
#

fileobj = Dataset(file0, mode='r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fill_value = fileobj.variables['sst']._FillValue.astype(float)
scale = fileobj.variables['sst'].scale_factor.astype(float)
offset = fileobj.variables['sst'].add_offset.astype(float)
fileobj.close()

#
# Size of block variables
#

years = np.unique(year)
NB = len(years)

#
# initialize some variables
#

alpha = 0.05
i_which = range(359,720+1,1) #,10)
j_which = range(150,360+1,1) #,10)
#i_which = range(0,X,4)
#j_which = range(0,Y,4)
lon_map = lon[i_which]
lat_map = lat[j_which]
X = len(lon_map)
Y = len(lat_map)
SST = {}
keys = ['TMM', 'TMX', 'TMN', 'TAM', 'TAX', 'TAN', 'TSpX', 'TSpN', 'TSuX', 'TSuN','TAuX', 'TAuN','TWiX', 'TWiN', 'TA90', 'TA10', 'TSp90', 'TSp10', 'TSu90', 'TSu10','TAu90', 'TAu10','TWi90', 'TWi10']
for key in keys:
    SST[key] = np.zeros((Y, X, NB))

SPR = np.in1d(month, np.array([10,11,12]))
SUM = np.in1d(month, np.array([1,2,3]))
AUT = np.in1d(month, np.array([4,5,6]))
WIN = np.in1d(month, np.array([7,8,9]))

#
# loop through locations
#

icnt = 0
for i in i_which:
    print icnt+1, 'of', len(lon_map)
#   load SST
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + str(i+1).zfill(4) + '.mat')
    sst_ts = matobj['sst_ts']
        #if np.logical_not(np.isfinite(sst_ts[j,:].sum())): # check for land
        #    jcnt += 1
        #    continue
#   Count number of MHWs of each length
        #mhws, clim = mhw.detect(t, sst_ts[j,:])
        #mhwBlock = mhw.blockAverage(t, mhws)
        # SST mean, trend-map, and time series
    for yr in range(len(years)):
        SST['TMM'][:,icnt,yr] = np.mean(sst_ts[j_which,:][:,year==years[yr]], axis=1)
        SST['TMX'][:,icnt,yr] = np.max(sst_ts[j_which,:][:,year==years[yr]], axis=1)
        SST['TMN'][:,icnt,yr] = np.min(sst_ts[j_which,:][:,year==years[yr]], axis=1)
        #SST['TAM'][:,icnt,yr] = 
        #SST['TAX'][:,icnt,yr] = 
        #SST['TAN'][:,icnt,yr] = 
        SST['TSpX'][:,icnt,yr] = np.max(sst_ts[j_which,:][:,(year==years[yr])*SPR], axis=1)
        SST['TSuX'][:,icnt,yr] = np.max(sst_ts[j_which,:][:,(year==years[yr])*SUM], axis=1)
        SST['TAuX'][:,icnt,yr] = np.max(sst_ts[j_which,:][:,(year==years[yr])*AUT], axis=1)
        SST['TWiX'][:,icnt,yr] = np.max(sst_ts[j_which,:][:,(year==years[yr])*WIN], axis=1)
        SST['TSpN'][:,icnt,yr] = np.min(sst_ts[j_which,:][:,(year==years[yr])*SPR], axis=1)
        SST['TSuN'][:,icnt,yr] = np.min(sst_ts[j_which,:][:,(year==years[yr])*SUM], axis=1)
        SST['TAuN'][:,icnt,yr] = np.min(sst_ts[j_which,:][:,(year==years[yr])*AUT], axis=1)
        SST['TWiN'][:,icnt,yr] = np.min(sst_ts[j_which,:][:,(year==years[yr])*WIN], axis=1)
        SST['TA90'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,year==years[yr]], 90, axis=1)
        SST['TA10'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,year==years[yr]], 10, axis=1)
        SST['TSp90'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*SPR], 90, axis=1)
        SST['TSu90'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*SUM], 90, axis=1)
        SST['TAu90'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*AUT], 90, axis=1)
        SST['TWi90'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*WIN], 90, axis=1)
        SST['TSp10'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*SPR], 10, axis=1)
        SST['TSu10'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*SUM], 10, axis=1)
        SST['TAu10'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*AUT], 10, axis=1)
        SST['TWi10'][:,icnt,yr] = np.percentile(sst_ts[j_which,:][:,(year==years[yr])*WIN], 10, axis=1)
    # Up counts
    icnt += 1
    # Save data so far
    outfile = '/home/ecoliver/Desktop/data/NESP/SST_annual_extremes_Aus'
    np.savez(outfile, lon_map=lon_map, lat_map=lat_map, years=years, SST=SST, i_which=i_which, j_which=j_which)



