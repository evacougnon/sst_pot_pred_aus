'''
    generate and save data to then use in plot_mapX.py
    data correspond to: number of days above the 90th percentile
    or 10th per pixels around Australia as well as the threshold 
    value for each pixel 

	firstly using the whole time series, then using the years 
	  corresponding to a specific climate mode (will be region
	  dependant!!!!)


    Also plot the map of the number of Marine Heat (Cold) Waves
    event using the marine extreme tool box developped by Eric 
    Oliver based on the definition of MHW by Hobbay et al. 2016

   Author: Eva A Cougnon
   Created: Feb 2017
   Last modification:

'''

# load required modules

import numpy as np
from scipy import io # load matlab file
from scipy import stats
from datetime import date
import marineHeatWaves as mhw
import time as time

# Load data
pathroot = '/home/ecoliver/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'

matobj_tmp = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.0001.mat')
lat_tmp = matobj_tmp['lat']
lon_tmp = matobj_tmp['lon']
res = 0.25 # resolution of the data -- quarter degree

# usefull numbers
MinYear = 1982
MaxYear =2016
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,31).toordinal()+1)

# define region of study 
lat_min = -55
lat_max = 10.25
lon_min = 90
lon_max = 180.25
# find the closest index from the lat/lon defined above
lat_min_id=np.nonzero((lat_tmp>(lat_min-res)) & (lat_tmp<(lat_min)))
lat_max_id=np.nonzero((lat_tmp>(lat_max-res)) & (lat_tmp<(lat_max)))
lon_min_id=np.nonzero((lon_tmp>(lon_min-res)) & (lon_tmp<(lon_min)))
lon_max_id=np.nonzero((lon_tmp>(lon_max-res)) & (lon_tmp<(lon_max)))

# allocate memory
lat_map = lat_tmp[lat_min_id[0][0]:lat_max_id[0][0]+1]
lat_id = range(lat_min_id[0][0],lat_max_id[0][0]+1,1)
lon_map = lon_tmp[lon_min_id[0][0]:lon_max_id[0][0]+1]
lon_id = range(lon_min_id[0][0],lon_max_id[0][0]+1,1)
SST_X = {}
keys = ['Np90','Np10','Tp90','Tp10','T_max','T_min','T_mean','T_std']
#keys = ['N_MHW','N_MCW','Dur_MHW','Dur_MCW']
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
Np90 -- number of days >= p90 threshold
Np10 -- number of days <= p10 threshold
N_MHW -- number of Marine Heat Waves events (following Hobbay's definition
N_MCW -- number of Marine Cold Waves events
Dur_MHW -- total duration of all the MHWs events during the studied period
Dur_MCW -- total duration of all the MCWs events during the studied period
'''
for key in keys:
    SST_X[key] = np.zeros((len(lat_map),len(lon_map)))

# basics stats
# going through the files
icnt_id = 0 # counting index
for i in lon_id:
#i = lon_id[0]
    t_lon = time.time()
    print(icnt_id+1, 'of', len(lon_id))
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' \
                        + str(i+1).zfill(4) + '.mat')
    sst_ts = matobj['sst_ts']

  # find p90 and p10 threshold as well as max and min
    SST_X['Tp90'][:,icnt_id] = np.percentile(sst_ts[lat_id,:],90,axis=1)
    SST_X['Tp10'][:,icnt_id] = np.percentile(sst_ts[lat_id,:],10,axis=1)
    SST_X['T_max'][:,icnt_id] = np.nanmax(sst_ts[lat_id,:],axis=1)
    SST_X['T_min'][:,icnt_id] = np.nanmin(sst_ts[lat_id,:],axis=1)  
    SST_X['T_mean'][:,icnt_id] = np.nanmean(sst_ts[lat_id,:],axis=1)
    SST_X['T_std'][:,icnt_id] = np.nanstd(sst_ts[lat_id,:],axis=1)
  
  # for each point count the number of times (days) the temperature
  # is above (or equal) p90 or below (or equal) p10 threshold
    jcnt_id = 0
    for j in lat_id:
#        t_lat = time.time()
#        print('lat',jcnt_id+1, 'of', len(lat_id))
       
        cnt_p90_tmp = np.zeros(len(sst_ts[0,:]))
        cnt_p10_tmp = np.zeros(len(sst_ts[0,:])) 
        cnt_p90_tmp[sst_ts[j,:] >= SST_X['Tp90'][jcnt_id,icnt_id]] = 1
        cnt_p10_tmp[sst_ts[j,:] <= SST_X['Tp10'][jcnt_id,icnt_id]] = 1
        SST_X['Np90'][jcnt_id,icnt_id] = np.sum(cnt_p90_tmp)
        SST_X['Np10'][jcnt_id,icnt_id] = np.sum(cnt_p10_tmp)
        '''
## marine heat/cold waves
        temp = sst_ts[j,:]
# checking if continent or not -- mhw does not handle empty matrix
        if ((np.isnan(np.min(temp)) == False) & (np.isnan(np.max(temp)) == False)):
# apply mhw definition
            mhws, clim_h = mhw.detect(dtime, temp, climatologyPeriod=[1982,2016], \
                                      coldSpells=False)
            SST_X['N_MHW'][jcnt_id,icnt_id] = mhws['n_events']
            SST_X['Dur_MHW'][jcnt_id,icnt_id] = np.sum(mhws['duration'])

            mcws, clim_c = mhw.detect(dtime, temp, climatologyPeriod=[1982,2016], \
                                  coldSpells=True)
            SST_X['N_MCW'][jcnt_id,icnt_id] = mcws['n_events']
            SST_X['Dur_MCW'][jcnt_id,icnt_id] = np.sum(mcws['duration'])
         '''  
#        elapsed_lat = time.time() - t_lat
#        print('elapsed time for each lat:', elapsed_lat)
        jcnt_id = jcnt_id + 1

    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt_id = icnt_id + 1

## save files
outfile = '/home/ecougnon/ana/SST_map_numberX_Aus'
#SST_map_mhw_Aus_good'
np.savez(outfile, lat_map=lat_map, lon_map=lon_map, SST_X=SST_X, \
         lat_id=lat_id, lon_id=lon_id)
print('saved')




