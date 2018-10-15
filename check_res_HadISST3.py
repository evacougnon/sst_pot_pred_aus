'''
    HadSST3 (5x5 monthly) to check when observationo were good enough for 
    good data in the region
     
    Plot the time series of the proportion of valuid data points in time compare
        to the total number of possible points in the area
    As well as the level of uncertainties?!
'''
# load required modules

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fname = '/home/ecougnon/Desktop/WorkingFigures/HadSST3/HadSST.3.1.1.0.measurement_and_sampling_uncertainty.nc'

# define indexes for lat lon of one pixel
lat_px_min = 11 # deg N
lat_px_max = -56
lon_px_min = 89 #150 # deg E
lon_px_max = 180

ds = xr.open_dataset(fname)['uncertainty']. \
        sel(latitude=slice(lat_px_max,lat_px_min), \
            longitude=slice(lon_px_min,lon_px_max)) #. \
#        stack(loc=('latitude','longitude'))
# total number of possible point on the grid 
# included land points! But as we'll look at a ratio
# that should not be a problem! we'll just never reach 1
tot_pts = len(ds.longitude)*len(ds.latitude)

ds_count = ds.count(dim=('latitude','longitude')) #'loc')
ds_mean = ds.mean(dim=('latitude','longitude'))
oceanMax_pts = ds_count.max()

ds_proportion = ds_count/oceanMax_pts

#ds_count.plot()
fig, axes = plt.subplots(2,1)
axes
ds_proportion.plot(ax=axes[0])
plt.grid()
plt.ylabel('proportion of point covered compared to the max')
ds_mean.plot(ax=axes[1])
plt.ylabel('mean uncertainties for the whole region')


plt.show(block=False)

