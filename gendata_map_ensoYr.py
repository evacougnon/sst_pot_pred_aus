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
   Last modification: May 2017 (update +ve, -ve limit of ENSO)

'''

# load required modules

import numpy as np
from scipy import io # load matlab file
from scipy import stats
from datetime import date
from calendar import monthrange
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
SST_enso = {}
keys = ['Tp90','Tp10','T_max','T_min']
#keys = ['N_MHW','N_MCW','Dur_MHW','Dur_MCW']
'''
Tp90 -- temperature threshold for p90
Tp10 -- temperature threshold for p10
N_MHW -- number of Marine Heat Waves events (following Hobbay's definition
N_MCW -- number of Marine Cold Waves events
Dur_MHW -- total duration of all the MHWs events during the studied period
Dur_MCW -- total duration of all the MCWs events during the studied period
'''
for key in keys:
    SST_enso[key] = np.zeros((len(lat_map),len(lon_map)))

## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 
file_enso = np.genfromtxt('/home/ecougnon/data/enso/MEI_index.txt', \
                        skip_header=10, skip_footer = 30, delimiter='\t')
str_id = np.nonzero((file_enso[:,0]>(MinYear-1)) \
                    & (file_enso[:,0]<(MinYear+1)))

enso_id = np.zeros(NumYears*12) # the idexes to use for our study (from 1982)
dtime_enso = [None] * (NumYears*12)
enso_daily = np.empty(NumDays)
k=0
l=0
d=0
for yy in np.arange(str_id[0][0],len(file_enso[:,0])):
    for mm in np.arange(1,12+1):
        enso_id[k] = file_enso[yy,mm]
        dtime_enso[k] = date(MinYear+l,mm,1).toordinal()
        enso_daily[d:d+monthrange(MinYear+l,mm)[1]] = enso_id[k]
        k = k + 1
        d = d + monthrange(MinYear+l,mm)[1]
    l = l + 1
enso_p_id = np.nonzero(enso_daily>=0.75)
enso_n_id = np.nonzero(enso_daily<=-0.75)
enso_neu_id = np.nonzero((enso_daily>-0.75) & (enso_daily<0.75))

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
# SSTs for positive ENSO
    #sst_enso = sst_ts
    sst_enso = sst_ts.copy()
    sst_enso[:,enso_p_id] = np.nan
    sst_enso[:,enso_neu_id] = np.nan

  # find p90 and p10 threshold as well as max and min
    SST_enso['Tp90'][:,icnt_id] = np.nanpercentile(sst_enso[lat_id,:],90,axis=1)
    SST_enso['Tp10'][:,icnt_id] = np.nanpercentile(sst_enso[lat_id,:],10,axis=1)
    SST_enso['T_max'][:,icnt_id] = np.nanmax(sst_enso[lat_id,:],axis=1)
    SST_enso['T_min'][:,icnt_id] = np.nanmin(sst_enso[lat_id,:],axis=1)  
  
  # for each point count the number of times (days) the temperature
  # is above (or equal) p90 or below (or equal) p10 threshold
    '''
    jcnt_id = 0
    for j in lat_id:
#        t_lat = time.time()
#        print('lat',jcnt_id+1, 'of', len(lat_id))
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
#           
#        elapsed_lat = time.time() - t_lat
#        print('elapsed time for each lat:', elapsed_lat)
        jcnt_id = jcnt_id + 1
        '''
    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt_id = icnt_id + 1

## save files
outfile = '/home/ecougnon/ana/SST_map_nina075Yr_Aus'
#SST_map_mhw_Aus_good'
np.savez(outfile, lat_map=lat_map, lon_map=lon_map, SST=SST_enso, \
         lat_id=lat_id, lon_id=lon_id)
print('saved')




