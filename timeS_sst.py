'''
 script to read the SST data for one pixel around Australia
  (also for a defined region?! -- average for the region)
  and save in sst_ts the SST time serie and plot it
  from the Daily 1982-2016 NOAA OI SST v2 "hi-res" 
   (https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html), 
    location on sverdrup: /home/ecoliver/data/sst/noaa_oi_v2/avhrr/

    Created: Fed 2017
    Author: EA Cougnon
    under development

    '''

# load modules
import numpy as np
from scipy import io
import netCDF4 as nc

#import matplotlib
#matplotlib.use("TkAgg") % used (once) to allow figures display
from matplotlib import pyplot as plt
# from calendar import monthrange # give number of days in a month
# but does not work for 1900...
import glob # finds files and directories
import datetime

# data folder name + saving folder
fname = '/home/ecoliver/data/sst/noaa_oi_v2/avhrr/'
sname = '/home/ecougnon/ana/'

# usefull numbers
MinYear = 1982
MaxYear =2016 
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
# define indexes for lat lon of the pixel or region
lat_min = 139 # -55.125N
lat_max = 400 #400 # 10.125N
lon_min = 359 #359 # 90E
lon_max = 720 # 180E

# allocate space for variables
sst_ts=np.empty(NumDays)  # time series SST set to nans
sst_ts.fill(np.nan)
dtime = np.empty(NumDays)
dtime.fill(np.nan) # NB: time is in days since 1978-01-01 00:00:00

# loop over all the years 
k=0
for year in range(MinYear,MaxYear+1):
    print(year)
# load data looping through all the daily files
    filenames = sorted(glob.glob(fname + str(year) + '/*.nc'))
    for f in filenames:
        data = nc.Dataset(f,'r')
# check if averaging in needed (for a region) or look at one point
        if lat_max > lat_min and lon_max > lon_min:
            sst_ts[k,] = np.squeeze(np.mean(data.variables['sst'] \
			            [:,:,lat_min:lat_max,lon_min:lon_max]))
        elif lat_max < lat_min or lon_max < lon_min:
            print('issue with lat/lon indexes definition')
        else:
            sst_ts[k,] = np.squeeze(data.variables['sst'][:,:,lat_min,lon_min])
        k = k+1

# create the matrix time useful for plotting
dtime = [datetime.date(MinYear,1,1) + datetime.timedelta(days=i) \
         for i in range(NumDays)]

# save the output file: sst_ts and dtime
np.savez(sname + 'ts_sst_oz',dtime=dtime,sst_ts=sst_ts)
# note: the file will be read as a np.lib with npzfile type do:
#       data = np.load('tmp_test.npz')
#       data.files
#   to read the name of the variables in the file and call them as:
#       data['dtime']


## plot time series
plot_ts = plt.plot(dtime,sst_ts)
plt.gcf().autofmt_xdate() # make the x axis more readable 
plt.show()


