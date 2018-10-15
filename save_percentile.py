'''
	calc and save percentile from big files

'''

import numpy as np
import xarray as xr

outfile = '/home/ecougnon/ana/OTEs_p10_p90_Aus'
fname = '/home/ecougnon/ana/SSTa_daily_Aus_shift_time_dim.nc'
lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
tim = xr.open_dataset(fname)['time']
SSTa = xr.open_dataset(fname)['SSTa']
SSTa = SSTa.sel(time=tim, lat=lat, lon=lon)
SSTa_p10 = np.nanpercentile(SSTa,10,axis=0)
SSTa_p90 = np.nanpercentile(SSTa,90,axis=0)

np.savez(outfile, SSTa_p10=SSTa_p10, SSTa_p90=SSTa_p90, lat=lat, \
         lon=lon, tim=tim)



