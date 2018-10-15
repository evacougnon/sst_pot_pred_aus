'''
 script to read the SST data for Australian region
 from the Daily 1982-2016 NOAA OI SST v2 "hi-res" 
 (https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html), 
 location on sverdrup: /home/ecoliver/data/sst/noaa_oi_v2/avhrr/

Created: Fed 2017
Author: EA Cougnon

'''

# load modules
import numpy as np
from scipy import io
import netCDF4 as nc

import matplotlib
#matplotlib.use("TkAgg") % used (once) to allow figures display
from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

# data folder name
fname = '/home/data/sst/noaa_oi_v2/avhrr/1982/'

# define indexes for lat lon Australia
lat_min = 145 #139 # -55.125N
lat_max = 215 #400 # 10.125N
lon_min = 359 # 90E
lon_max = 720 # 180E

# load data 
data = nc.Dataset(fname + 'avhrr-only-v2.19821120.nc','r')
lat = data.variables['lat'][lat_min:lat_max]
lon = data.variables['lon'][lon_min:lon_max]
sst = np.squeeze(data.variables['sst'][:,:,lat_min:lat_max,lon_min:lon_max])
# sst_anom = data.variables['anom'][:] # daily sst anomalies

# reshape lat/lon as matrix
[lat_,lon_] = sst.shape
lat_tmp = lat[:,np.newaxis]
lats = np.tile(lat_tmp,(1,lon_))
lons = np.tile(lon,(lat_,1))


# plot sst
plt.figure(figsize=(8,8)) #10,4))
plt.clf()
domain = [-50, 140, -35, 160] #[-55, 90, 10, 180]
domain_draw = [-50, 140, -35, 160] #[-55, 90, 10, 180]
dlon = 20
dlat = 10
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.drawcoastlines(zorder=1)
proj.drawcountries(zorder=1)
proj.drawparallels(np.arange(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
proj.drawmeridians(np.arange(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
lonproj, latproj = proj(lons, lats)
plt.contourf(lonproj, latproj, sst[:,:], levels=np.arange(8,19,1), cmap=plt.cm.viridis)
h=plt.colorbar()
plt.show()







