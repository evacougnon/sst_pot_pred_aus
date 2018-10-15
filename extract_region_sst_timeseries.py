'''
    extract from the global daily SSTs files
    the sst for the wanted region (define lat/lon 
    min/max and do the average of the domain of 
    interest by doing a simple spatial average
    or an area(weighted) average -- see the "WHAT"
    option.

    Created: Mar 2018
    Author: Eva C.
    Last Modif:
'''


# load required modules

import numpy as np
import xarray as xr
import time as time
from geopy.distance import vincenty # distance lat/lon

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo
import eac_useful as eac


# define indexes for lat lon of one pixel
lat_px_min = -45 #-45 #-46 # deg N
lat_px_max = -37 #-37 #-26
lon_px_min = 147 #139 #147 #135 #150 # deg E
lon_px_max = 155 #147 #155 #174

header = '/home/data/sst/noaa_oi_v2/new_avhrr/'

# Choose whether it is area-averaged (A)
# or just spatially averaged (S) considering a regular
# lat/lon grid ( degree for instance)
# A_eo is following Eric's method to calculate lat/lon distance
WHAT = 'S' # 'S'

if WHAT == 'S': # spatial average
    outfile = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/TasSea_Oliver2017_bdr.nc'
# for the TAS Sea box (with 150oE western bdr) it takes: less than 20 minutes 
    t_read = time.time()
    res = 0.25
    ds = xr.open_mfdataset(header + 'sst.day.mean.????.v2.nc')['sst']. \
            sel(lat=slice(lat_px_min-res,lat_px_max+res), \
                lon=slice(lon_px_min-res,lon_px_max+res), \
                time=slice('1982-01-01','2018-05-19')).mean(dim=('lat','lon'))
    ds_oisst_d = xr.concat(ds, 'time').squeeze()
    ds_oisst_d.to_netcdf(outfile)
    elapsed_1 = time.time() - t_read
    print('elapsed time for reading/concat/save the region and doing the spatial mean: ', \
          elapsed_1)

elif WHAT == 'A': # area-average
    outfile = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/TasSea_Oliver2017_AreaAve_updated_May.nc'
#TasSea_135_AreaAve_masked_updated.nc'
    res = 0.25 # resolution of the lat/lon grid (in degree)
    t_read = time.time()
# calc the dx/dy of each cells, take one extra point along the 
# boundary to be able to calculate the area of the boundary cells
# open_dataset does not load the file but only keeps the link in memory
    ds_tmp = xr.open_dataset(header + 'sst.day.mean.1982.v2.nc')['sst']. \
                sel(lat=slice(lat_px_min-res,lat_px_max+res), \
                    lon=slice(lon_px_min-res,lon_px_max+res), \
                    time=slice('1982-01-01','1982-12-31'))
# load lat/lon
    lon = ds_tmp.lon
    lat = ds_tmp.lat
# get the distance between each grid points in meters
# dx abd dy are of dimension (lat,lon)
    dx, dy = eac.dxdy(lat,lon)

# load in memory the link to all the files of interest
    ds = xr.open_mfdataset(header + 'sst.day.mean.????.v2.nc')['sst']. \
            sel(lat=slice(lat_px_min,lat_px_max), \
                lon=slice(lon_px_min,lon_px_max), \
                time=slice('1982-01-01','2018-05-19'))
# The open_mfdataset (opens multiple files) does not accept a masked array 
# as the single file command (open_dataset) does!! Needs a matrix of 
# zeros and ones, or nans and ones!
# create a mask array to NaN the land
    mask_land = np.ones((len(ds.lat),len(ds.lon)))
    tmp = np.squeeze(ds[0,:,:])
    mask_land = np.ma.masked_where(np.isnan(tmp)==True, mask_land)
# calc area by masking (ignoring) land points
    area = dx * dy
    area_mask = area.copy()    
    area_mask[np.where(mask_land==False)] = np.nan
# divide each ocean-cell area by the total ocean area 
    weight = area_mask/np.nansum(area_mask)
    tot_weight = np.nansum(weight) # quick check, needs to be 1
    '''
    # what I had before!!!
    # confusion between the True/False, I used to work in matlab
    # with 0 and 1 instead of True and False, somehow I had in mind
    # that 1 (ocean in that case) was True! However for a masked
    # array in python True means that it is the mask (so 0 in by matlab
    # way of thinking!!!...)
    area = dx * dy
    area[np.where(mask_land==True)] = np.nan
    area_tot = np.nansum(area)
    weight = area/area_tot
    tot_weight = np.nansum(weight) # This total weigh check wasn't
    			# equals to 1! I don't know how I didn't 
			# pick this error! I was probably checking 
			# too often the open single file command to 
			#test my code!

    '''
# orignal ds is of dimension (time, lat, lon), 
# ds*weight automatically multiply on the lat,lon dimension
# checked several time steps using mfdataset on 1982 and 198? files
# Further note: this below line could also be part of a 
# 'prepross' argument in the mfdataset command
# for more efficiency -- think about implementing it
# in the future, can be very useful for wider area
# check the manual of the command
    ds = ((ds*weight).sum(dim=('lon','lat'))) / tot_weight
    '''    
    # testing script!!!!!
    ds_test = np.array(ds)
    ds_out = np.NaN*np.zeros((len(ds.time)))
    for tt in np.arange(0,len(ds.time)):
        tmp = np.NaN*np.zeros((len(ds.lat),len(ds.lon)))
        tmp[:,:] = ds_test[tt,:,:]*weight[:,:]
        ds_out[tt] = np.nansum(tmp)
   
    '''
    ds_oisst_d = xr.concat(ds, 'time').squeeze()
    ds_oisst_d.to_netcdf(outfile)
    
    elapsed_1 = time.time() - t_read
    print('elapsed time for reading/concat/save the region and doing the spatial mean: ', \
          elapsed_1)

elif WHAT == 'A_eo':
    print('needs update!!!! -- wrong tot_weight!!! see above!')
    '''
    outfile = '/home/ecougnon/ana/MHW_paper_SarahKirkpatrick/TasSea_135_AreaAve_masked_eo.nc'
    res = 0.25 # resolution of the lat/lon grid (in degree)
    t_read = time.time()
    ds_tmp = xr.open_dataset(header + 'sst.day.mean.1982.v2.nc')['sst']. \
                sel(lat=slice(lat_px_min-res,lat_px_max), \
                    lon=slice(lon_px_min-res,lon_px_max), \
                    time=slice('1982-01-01','1982-12-31'))
# calc distance following Eric's method
    lon = ds_tmp.lon
    lat = ds_tmp.lat
    dx, dy = eo.dxdy(lon,lat) # in meters
# dask array for all the files
    ds = xr.open_mfdataset(header + 'sst.day.mean.????.v2.nc')['sst']. \
            sel(lat=slice(lat_px_min,lat_px_max), \
                lon=slice(lon_px_min,lon_px_max), \
                time=slice('1982-01-01','2018-02-28'))
# create a ask to NaN the land
    mask_land = np.ones((len(ds.lat),len(ds.lon)))
    tmp = np.squeeze(ds[0,:,:])
    mask_land = np.ma.masked_where(np.isnan(tmp)==True, mask_land)
# calc area by masking (ignoring) land points
    area = (dx*1000)*(dy*1000) # get meters
    area_tot = area.copy()    
    area_tot[np.where(mask_land==False)] = np.nan
    area_tot = np.nansum(area_tot)
    weight = area/area_tot
    tot_weight =np.nansum(weight)

    ds = ((ds*weight).sum(dim=('lon','lat'))) / tot_weight

    ds_oisst_d = xr.concat(ds, 'time').squeeze()
    ds_oisst_d.to_netcdf(outfile)
    elapsed_1 = time.time() - t_read
    print('elapsed time for reading/concat/save the region and doing the spatial mean: ', \
          elapsed_1)
    '''


print('saved in', outfile)


