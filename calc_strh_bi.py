'''
   read mslp data downloaded from 
   https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2.monolevel.mm.html

   and calc the subtropical ridge tasman high index following Andrew Marshall 
   technique (Marshall et al 2013 DOI 10.1007/s00382-013-2016-1)

   lat and lon limits are defined followin his paper:
   "STRH index, [...] MSLP data averaged over the box region 
   (150E–165E, y?S ± 7.5?). Here, yS is the climatological central latitude 
   of the subtropical ridge, determined for each month over our analysis period, 
   which varies in its annual cycle from 27.5S in August to 40Sin February over 
   the Tasman Sea.

   CR20 reanalysis goes only until Decembre 2012!

'''

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

sname_strh = '/home/ecougnon/strh_index'
sname_bi_140 = '/home/ecougnon/blocking_index_140'
sname_bi_160 = '/home/ecougnon/blocking_index_160'

# load data
fname_strh = '/home/ecougnon/prmsl.mon.mean.nc'
lat = xr.open_dataset(fname_strh)['lat']
lat = lat.sel(lat=slice(-20.,-50.))
lon = xr.open_dataset(fname_strh)['lon']
lon = lon.sel(lon=slice(150.,165.))
time_ts = xr.open_dataset(fname_strh)['time']
time_ts = time_ts.sel(time=slice('1982-01-01','2012-12-01'))
mslp_ini = xr.open_dataset(fname_strh)['prmsl']
# to use in STRH calc
mslp_tmp = mslp_ini.sel(time=slice('1982-01-01','2012-12-01'), lat=slice(-10.,-60.),\
                         lon=slice(150.,165.))
# to use to identify the central ridge
mslp_ = mslp_ini.sel(time=slice('1982-01-01','2012-12-01'), lat=slice(-20.,-50.), \
                                lon=slice(150.,165.))

# calc the monthly climatology of the central latitude of the 
# subtropical ridge
zonal_mslp = np.mean(mslp_,axis=2)

# monthly latitude of the max mslp to define the  central ridge
lat_max_ts = np.zeros(len(time_ts))
for tt in range(0,len(time_ts)):
    max_id = np.where(zonal_mslp[tt,:] == np.max(zonal_mslp[tt,:]))
    lat_max_ts[tt] = lat[max_id[0][0]] # take the first index when several max_id
# calc the monthly climatology of the latitude
lat_max_clim = np.zeros(12)
for tt in range(0,12):
    lat_max_clim[tt] = np.mean(lat_max_ts[range(tt,len(time_ts),12)])

# calc STRH index in Pa -- not anomalies
strh = np.empty((len(time_ts)))
strh.fill(np.nan)
k=0
for tt in range(0,12):
    lat_tmp = lat_max_clim[tt]
    tmp = mslp_tmp.sel(lat=slice(lat_tmp+7.5,lat_tmp-7.5))
    strh[range(tt,len(time_ts),12)]=np.mean(np.mean(tmp[range(tt,len(time_ts),12), \
                                                        :,:], axis=1),axis=1)

np.savez(sname_strh, strh = strh, time_ts=time_ts)


# blocking index
fname_bi = '/home/ecougnon/uwnd.mon.mean.nc'
time_bi = xr.open_dataset(fname_bi)['time']
time_bi = time_bi.sel(time=slice('1982-01-01','2012-12-01'))
uwnd = xr.open_dataset(fname_bi)['uwnd']

uwnd_140 = uwnd.sel(time=slice('1982-01-01','2012-12-01'), level=500.0, lon=140.)
u_25 = uwnd_140.sel(lat=-25., method='nearest')
u_30 = uwnd_140.sel(lat=-30., method='nearest')
u_40 = uwnd_140.sel(lat=-40., method='nearest')
u_45 = uwnd_140.sel(lat=-45., method='nearest')
u_50 = uwnd_140.sel(lat=-50., method='nearest')
u_55 = uwnd_140.sel(lat=-55., method='nearest')
u_60 = uwnd_140.sel(lat=-60., method='nearest')

BI_140 = 0.5*(u_25 + u_30 - u_40 - 2*u_45 - u_50 + u_55 + u_60)

uwnd_160 = uwnd.sel(time=slice('1982-01-01','2012-12-01'), level=500.0, lon=160.)
u_25_ = uwnd_160.sel(lat=-25., method='nearest')
u_30_ = uwnd_160.sel(lat=-30., method='nearest')
u_40_ = uwnd_160.sel(lat=-40., method='nearest')
u_45_ = uwnd_160.sel(lat=-45., method='nearest')
u_50_ = uwnd_160.sel(lat=-50., method='nearest')
u_55_ = uwnd_160.sel(lat=-55., method='nearest')
u_60_ = uwnd_160.sel(lat=-60., method='nearest')

BI_160 = 0.5*(u_25_ + u_30_ - u_40_ - 2*u_45_ - u_50_ + u_55_ + u_60_)

np.savez(sname_bi_140, BI_140 = BI_140, time_bi=time_bi)
np.savez(sname_bi_160, BI_160 = BI_160, time_bi=time_bi)





