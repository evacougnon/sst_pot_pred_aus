'''
	generate the Ningaloo Nino Index from the monthly
	SSTa from the NOAA OISST data set. Anomalie calculated 
	from the 1982-2016 time series

	Author: Eva C
	Created: Nov 2017
'''
# libriries
import numpy as np
import xarray as xr
import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_useful as eac
from matplotlib import pyplot as plt
import matplotlib

figfile = '/home/ecougnon/ana/Ningaloo_Nino_index.png'
outfile = '/home/ecougnon/data/NingalooNino_index'


# load data
fname = '/home/ecougnon/ana/SSTa_monthly_extremes_Aus.nc'
lat = xr.open_dataset(fname)['lat']
lat = lat.sel(lat=slice(-32,-22))
lon = xr.open_dataset(fname)['lon']
lon = lon.sel(lon=slice(110,120))
tim = xr.open_dataset(fname)['time']
# for HadISST
#tim = tim.sel(time=slice('1871-01-01','2016-01-01'))
SST = xr.open_dataset(fname)
SST = SST['TMM'].sel(lat=lat, lon=lon)

NingN_Ind = np.squeeze(np.apply_over_axes(np.nanmean, SST, (2,1)))

# saving
np.savez(outfile, NingN_Ind=NingN_Ind, tim=tim, SST=SST, lat=lat, lon=lon)

# plotting time series
plt.figure(figsize=(13,7))
plt.plot(tim, NingN_Ind)
plt.plot(tim[(12*4/2):-(12*4/2)+1], eac.moving_average(NingN_Ind, 12*4))
plt.grid()
plt.title('Ningaloo Nino index from the monthly SSTa -- 22-32S and 110-coast')
plt.ylabel('SSTa anomalies')
plt.legend(['SSTa', '4 year moving average'])

plt.savefig(figfile, bbox_inches='tight',dpi=300)
plt.show()



'''
# plotting -- maps
domain = [-35, 100, -20, 120] #[-55, 90, 10, 180] 
domain_draw = [-35, 100, -20, 120] #[[-50, 90, 10, 180] 
dlat = 2 #10
dlon = 2 #30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure()
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, np.nanmean(SST, axis=0), levels=np.arange(-5,5+1,1), \
             cmap=plt.cm.YlOrBr)
'''











