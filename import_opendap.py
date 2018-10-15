'''
   import data with opendap
   in this case MSLP from CR20
'''

import netCDF4 
import xarray as xr
import numpy as np

url = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/20thC_ReanV2/Monthlies/monolevel/prmsl.mon.mean.nc?lat[0:1:90],lon[0:1:179],time_bnds[0:1:1703][0:1:1],time[0:1:1703],prmsl[0:1:1703][0:1:90][0:1:179]'

dataset = netCDF4.Dataset(url)
mslp = dataset.variables['prmsl']
lat = dataset.variables['lat']
lon = dataset.variables['lon']
tim = dataset.variables['time']

#t_file = netCDF4.Dataset('test.nc','w',format='NETCDF4')


#from pydap.client import open_url
#dataset2 = open_url('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/20thC_ReanV2/Monthlies/monolevel/prmsl.mon.mean.nc?lat[0:1:90],lon[0:1:179],time_bnds[0:1:1703][0:1:1],time[0:1:1703],prmsl[0:1:1703][0:1:90][0:1:179]', encoding='utf-8', errors='ignore')


MSLP = xr.Dataset({'mslp':(('time','lat','lon'),test)}, \
                  {'time': tim, 'lat': lat, 'lon':lon})



