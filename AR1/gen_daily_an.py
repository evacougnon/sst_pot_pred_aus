'''
    do the anomalies on the daily by removing the seaon
    using xarray 

    Created: Oct 2018
    Author: Eva C.
'''

import numpy as np
import pandas as pd
import xarray as xr
import dask.array
import time as time

import sys
sys.path.insert(0,'../libraries/')
import eric_oliver as eo

# saving file -- with nc, does not allow to overwrite!!!
outfile = '/data/home/evac/NESP/PotPred/SST_daily_trend_Aus.nc'
outfile2 = '/data/home/evac/NESP/PotPred/SSTa_daily_trend_Aus.nc'

header = '/data/home/oliver/data/sst/noaa_oi_v2/avhrr/'
lat_min = -55 
lat_max = 10
lat_bin = 1 
lon_min = 90 
lon_max = 180 
lon_bin = 1 

tim = pd.date_range('1982-01-01','2016-12-31')
# remove the last day of the year when leap year
# Need a function!!!....
tim = tim[tim !='1984-12-31']
tim = tim[tim !='1988-12-31']
tim = tim[tim !='1992-12-31']
tim = tim[tim !='1996-12-31']
tim = tim[tim !='2000-12-31']
tim = tim[tim !='2004-12-31']
tim = tim[tim !='2008-12-31']
tim = tim[tim !='2012-12-31']
tim = tim[tim !='2016-12-31']
'''
ds = xr.open_mfdataset(header + 'sst.day.mean.????.nc')['sst']. \
        sel(time=tim, lat=slice(lat_min,lat_max), \
            lon=slice(lon_min,lon_max))
ds = xr.concat(ds, 'time').squeeze()
ds.to_netcdf(outfile)
'''

t_ini = time.time()

ds = xr.open_dataset(outfile)['sst']

SSTa = xr.Dataset({'SSTa':(('time','lat','lon'), np.zeros((ds.shape[0],ds.shape[1],ds.shape[2])))},coords = {'time':ds.time, 'lat': ds.lat.coords['lat'], 'lon':ds.lon.coords['lon']})
tot_pts = len(ds.lat)*len(ds.lon)
k = 0
for ll in range(0,len(ds.lat)):
    for ln in range(0,len(ds.lon)):
#        t_key = time.time()
        tmp_sea = eo.deseason_harmonic(np.array(ds[:,ll,ln]),4,365)
        SSTa['SSTa'][:,ll,ln] = np.array(tmp_sea)[:,0]
    tmp = (k+len(ds.lon))/tot_pts
    k = k+len(ds.lon)
    print('percentage of points processess:', tmp)
#        elapsed_key = time.time() - t_key
#        print('elapsed time for each key:', elapsed_key)

SSTa.to_netcdf(outfile2)

elapsed_all = time.time() - t_ini
print('elapsed time since start:', elapsed_all)



