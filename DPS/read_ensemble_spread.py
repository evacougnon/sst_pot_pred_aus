'''
    read and do some analyses with the ensembe spread
    to identify where the model gets SST extremes

'''


# load required modules

import numpy as np
import xarray as xr
import pandas as pd

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

# load data
header_out = '/home/ecougnon/data/DPS/forecast/plots/' #reanalysis/ETKF/'
header = '/home/ecougnon/data/DPS/forecast/' #reanalysis/ETKF/'
fname = header + 'temp_spread_enkf13.nc' #'SST_ensemble_spread_2002-2017.nc'

# define the region of interest
# Australian region in the spread file (not in lat/lon!)
lat_min = 45
lat_max = 160
lon_min = 1 #10
lon_max = 210 #110
yaxis_1 = xr.open_dataset(fname)['yaxis_1']
yaxis_1 = yaxis_1.sel(yaxis_1=slice(lat_min,lat_max))
xaxis_1 = xr.open_dataset(fname)['xaxis_1'] 
xaxis_1 = xaxis_1.sel(xaxis_1=slice(lon_min,lon_max))

ID = 66 
# ID=13 with the spread from enkf-9 and starting in 2003
# ID = 14 with the spread from enkf-13 starting in 2003
# ID = 66 with the spread from enkf-13 starting in 2007

T = len(xr.open_dataset(fname)['Time'])
'''
trying to organise tim_vec readin ncdump -h fname (contains the information
with the names of the concatenated files

tim_vec = np.zeros(T)
tim_vec_1 = pd.date_range('2002-01-31','2008-01-31',name='time',freq='30D')
tim_vec_2

tim_vec[:len(tim_vec_1)] = tim_vec_1
tim_vec[len(tim_vec_1)] = pd.date('2008-02-27')
time_vec[len(tim_vec_1)+1:??? = tim_vec_2
'''
tim_vec = pd.date_range('2002-01-01','2010-10-31',name='time',freq='28D')
#tim_vec = pd.date_range('2002-01-01','2017-12-31',name='time',freq='29D')
#tim_vec = pd.date_range('2003-01-26','2017-12-31',name='time',freq='30D')


sst_spread = xr.open_dataset(fname)['temp'].sel(Time=slice(ID,T),xaxis_1=xaxis_1, \
                                                yaxis_1=yaxis_1, \
                                                zaxis_1=0).squeeze()
 #['sstb_an'] #['temp']


## plot time series of the spread
'''
74,74: point in the EAC -- near the coast ~ Sydney
'''
# define the number of points wanted:
N=4
# define the name of the location/point
keys = ['Sydney','NE_NZ','W_Aus','EqW_Pac']
test_point = np.zeros((N,len(tim_vec[ID:])))

test_point[0,:] = xr.open_dataset(fname)['temp'].sel(Time=slice(ID,T), \
                                                     xaxis_1=74, \
                                                     yaxis_1=74, \
                                                     zaxis_1=0).squeeze()
test_point[1,:] = xr.open_dataset(fname)['temp'].sel(Time=slice(ID,T), \
                                                     xaxis_1=105, \
                                                     yaxis_1=78, \
                                                     zaxis_1=0).squeeze()
test_point[2,:] = xr.open_dataset(fname)['temp'].sel(Time=slice(ID,T), \
                                                     xaxis_1=32, \
                                                     yaxis_1=77, \
                                                     zaxis_1=0).squeeze()
test_point[3,:] = xr.open_dataset(fname)['temp'].sel(Time=slice(ID,T), \
                                                     xaxis_1=105, 
                                                     yaxis_1=137, \
                                                     zaxis_1=0).squeeze()
'''
plt.figure(figsize=(10,4))
plt.plot(tim_vec[ID:], test_point[0,:])
plt.plot(tim_vec[ID:], test_point[1,:])
plt.plot(tim_vec[ID:], test_point[2,:])
plt.plot(tim_vec[ID:], test_point[3,:])
plt.legend(['Sydney', 'NE NZ','NW Aus','Eq W Pac'], fontsize=10)
plt.title('ensemble spread at specific points (enkf-13)', fontsize=14)
plt.grid()
plt.savefig(header_out +'ensemble_spread_enks13_points.png', \
            bbox_inches='tight', format='png', dpi=300)

#plt.show()
'''


sst_spread_mean = np.array(sst_spread.mean(dim='Time'))
sst_spread_mean[sst_spread_mean==0] = np.nan

sst_spread_max = np.array(sst_spread.max(dim='Time'))
sst_spread_max[sst_spread_max==0] = np.nan

## plot the mean
plt.figure(figsize=(11,9))
plt.clf()
ax = plt.subplot(211)
plt.contourf(xaxis_1, yaxis_1, sst_spread_mean, \
             levels=np.arange(0.2,0.9+0.05,0.05) , \
             cmap=plt.cm.viridis)
cb = plt.colorbar(ticks=np.arange(0.2,0.9+0.1,0.1),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.title('mean of the ensemble spread 2007-2010', fontsize=14)

ax = plt.subplot(212)
plt.contourf(xaxis_1, yaxis_1, sst_spread_max, \
             levels=np.arange(0.2,1.5+0.1,0.1) , \
             cmap=plt.cm.viridis)
cb = plt.colorbar(ticks=np.arange(0.2,1.5+0.1,0.2),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.title('max of the ensemble spread 2007-2010', fontsize=14)

#plt.savefig(header_out +'map_ensemble_spread_enkf13_mean_max_Pac.png', \
#            bbox_inches='tight', format='png', dpi=300)
plt.show()









