'''
    time series analyses for specific points chosen
    around Australia to compare to different climate
    mode -- starting with ENSO

   script based in some extend on ts_tmp.py script (testing script)


    Author: Eva A Cougnon
    Created: Mar 2017
    Last modified: May 2017 

'''


# load required modules

import numpy as np
import datetime 
from calendar import monthrange

from matplotlib import pyplot as plt
from matplotlib import mlab as mlab # Numerical python 
        # functions written for compatability with MATLAB 
        # commands with the same names
# usefull numbers
MinYear = 1982
MaxYear =2016
NumYears = MaxYear-MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
# time vector for plotting
dtime = [datetime.date(MinYear,1,1) + datetime.timedelta(days=i) \
         for i in range(NumDays)]
# define indexes for lat lon of one pixel
lat_px =  -5 #-30 #-42 #-30 # deg N
lon_px = 180 #112 #151 #112 # deg E

# fig files saving names
figfile_ts = '/home/ecougnon/ana/ts_ssta_enso_5-180.eps'
figfile_hist = '/home/ecougnon/ana/hist_ssta_enso_5-180.eps'

## whether or not using the SSTa_daily_Aus.nc
# check = 0 -- pick a file from the original SST.mat files
# check = 1 -- read a nc file 
check = 1

if (check ==0):  
    from scipy import io # load matlab file
# Load data
    pathroot = '/home/ecoliver/'
    header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
    matobj_tmp = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.0001.mat')
    lat_tmp = matobj_tmp['lat']
    lon_tmp = matobj_tmp['lon']
    res = 0.25 # resolution of the data -- quarter degree
# find the closest index from the lat/lon defined above
    lat_px_id=np.nonzero((lat_tmp>(lat_px-res)) & (lat_tmp<(lat_px)))
    lon_px_id=np.nonzero((lon_tmp>(lon_px-res)) & (lon_tmp<(lon_px)))
# load time series from the given point
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' \
                        + str(lon_px_id[0][0]).zfill(4) + '.mat')
    sst_ts = matobj['sst_ts'][lat_px_id[0][0],:]
elif (check == 1):
    import pandas as pd
    import xarray as xr
    fname = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
    lat = xr.open_dataset(fname)['lat']
    lat = lat.sel(lat=lat_px, method='nearest') 
    lon = xr.open_dataset(fname)['lon']
    lon = lon.sel(lon=lon_px, method='nearest') 
    tim = xr.open_dataset(fname)['time']
    sst_ts = xr.open_dataset(fname)['SSTa']
    sst_ts = sst_ts.sel(lat=lat,lon=lon)
## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 
file_enso = np.genfromtxt('/home/data/index/enso/MEI_index.txt', \
                        skip_header=10, skip_footer = 30, delimiter='\t')
# starting index along the MEI time series that starts in 1950, while
# our time series starts in 1982
str_id = np.nonzero((file_enso[:,0]>(MinYear-1)) \
                    & (file_enso[:,0]<(MinYear+1)))
enso = np.zeros(NumYears*12) # ENSO index values for our study (from 1982)
dtime_enso = [None] * (NumYears*12) # corresponding time vector for plotting
enso_daily = np.empty(NumDays)
k=0
l=0
d=0
for yy in np.arange(str_id[0][0],len(file_enso[:,0])):
    for mm in np.arange(1,12+1):
        enso[k] = file_enso[yy,mm]
        dtime_enso[k] = datetime.date(MinYear+l,mm,1)
        enso_daily[d:d+monthrange(MinYear+l,mm)[1]] = enso[k]
        k = k + 1
        d = d + monthrange(MinYear+l,mm)[1]
    l = l + 1
# save indexes where enso is +ve (p) and -ve (n)
enso_p_id = np.nonzero(enso_daily>=0.75)
enso_n_id = np.nonzero(enso_daily<=-0.75)
enso_neu_id = np.nonzero((enso_daily>-0.75) & (enso_daily<0.75))
# calc stats from the point time serie
sst_ts_mean = np.mean(sst_ts)
sst_ts_std = np.std(sst_ts)
sst_ts_p10 = np.percentile(sst_ts,10)
sst_ts_p90 = np.percentile(sst_ts,90)
sst_ts_med = np.median(sst_ts)
sst_ts_max = np.max(sst_ts)
sst_ts_min = np.min(sst_ts)

# calc stats from the point time serie
# considering +ve or -ve ENSO
sst_enso_n = sst_ts.copy() # to use only on numpy array
sst_enso_p = sst_ts.copy()
sst_enso_n[enso_p_id[0][:]] = np.nan
sst_enso_n[enso_neu_id[0][:]] = np.nan
sst_enso_p[enso_n_id[0][:]] = np.nan
sst_enso_p[enso_neu_id[0][:]] = np.nan

sst_enso_p_mean = np.nanmean(sst_enso_p)
sst_enso_p_std = np.nanstd(sst_enso_p)
sst_enso_p_p10 = np.nanpercentile(sst_enso_p,10)
sst_enso_p_p90 = np.nanpercentile(sst_enso_p,90)
sst_enso_p_med = np.nanmedian(sst_enso_p)
sst_enso_p_max = np.nanmax(sst_enso_p)
sst_enso_p_min = np.nanmin(sst_enso_p)

sst_enso_n_mean = np.nanmean(sst_enso_n)
sst_enso_n_std = np.nanstd(sst_enso_n)
sst_enso_n_p10 = np.nanpercentile(sst_enso_n,10)
sst_enso_n_p90 = np.nanpercentile(sst_enso_n,90)
sst_enso_n_med = np.nanmedian(sst_enso_n)
sst_enso_n_max = np.nanmax(sst_enso_n)
sst_enso_n_min = np.nanmin(sst_enso_n)


## PLOTTING
plt.figure(figsize=(13,19))
plt.clf()
# plot p10 threshold
ax1=plt.subplot(2,1,1)
plt.plot(dtime_enso,enso, color='k')
ax1.set_xlim([dtime_enso[0], dtime_enso[-1]])
ax1.set_ylim([-2, 3])
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
ax1.fill_between(dtime_enso, enso, where=enso >= 0.75, facecolor='red', alpha=0.5, \
                 interpolate=True)
ax1.fill_between(dtime_enso, enso, where=enso <= -0.75, facecolor='blue', alpha=0.5, \
                 interpolate=True)
plt.title('Multivariate ENSO Index', fontsize=16, y=1.02)
plt.grid()
#plt.gcf().autofmt_xdate() # make the x axis more readable 
ax2=plt.subplot(2,1,2)
plt.plot(dtime,sst_ts, color='0.7')
plt.plot(dtime,sst_enso_n, color='b', linestyle = ':')
plt.plot(dtime,sst_enso_p, color='r')
plt.axhline(y=sst_ts_p10, color='k', linestyle=':', linewidth=1, label='time series')
plt.axhline(y=sst_enso_n_p10, color='b', linestyle='-.', linewidth=1, label='-ENSO')
plt.axhline(y=sst_enso_p_p10, color='r', linestyle='--', linewidth=1, label='+ENSO')
plt.axhline(y=sst_ts_p90, color='k', linestyle=':', linewidth=1)
plt.axhline(y=sst_enso_n_p90, color='b', linestyle='-.', linewidth=1)
plt.axhline(y=sst_enso_p_p90, color='r', linestyle='--', linewidth=1)
ax2.legend()
plt.title('SSTa at 5S 180E -- East PNG', fontsize=16, y=1.02)
ax2.set_xlim([dtime[0], dtime[-1]])
ax2.set_ylim([-4, 4])
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
plt.grid()
plt.savefig(figfile_ts,bbox_inches='tight', format='eps', dpi=300)

#'''
# histogram
num_bins = 100 # in the future make it every 0.1 deg C
fig, ax = plt.subplots()
# the histogram of the data
#n, bins, patches = ax.hist([sst_enso_n[~np.isnan(sst_enso_n)]], num_bins, \
#normed=1, histtype='bar',stacked=True)
plt.hist(sst_enso_n[~np.isnan(sst_enso_n)], num_bins, normed = True, \
         color='b', alpha=0.5, stacked=True)
plt.hist(sst_enso_p[~np.isnan(sst_enso_p)], num_bins, normed = True, \
         color='r', alpha=0.5, stacked=True)
plt.hist(sst_ts, num_bins, normed = True, color='k', alpha=0.7, stacked=True)
# add a 'best fit' line
#y = mlab.normpdf(bins, sst_ts_mean, sst_ts_std)
#ax.plot(bins, y, '--')
ax.axvline(x=sst_ts_p10, ymin=0, ymax=1, linewidth=1, color='k', linestyle=':', label='time series')
ax.axvline(x=sst_ts_p90, ymin=0, ymax=1, linewidth=1, color='k', linestyle=':')
ax.axvline(x=sst_ts_mean, ymin=0, ymax=1, linewidth=1, color='k', linestyle=':')
ax.axvline(x=sst_enso_p_p10, ymin=0, ymax=1, linewidth=1, color='r', linestyle='--', label='+ENSO')
ax.axvline(x=sst_enso_p_p90, ymin=0, ymax=1, linewidth=1, color='r', linestyle='--')
ax.axvline(x=sst_enso_p_mean, ymin=0, ymax=1, linewidth=1, color='r', linestyle='--')
ax.axvline(x=sst_enso_n_p10, ymin=0, ymax=1, linewidth=1, color='b', linestyle='-.', label='-ENSO')
ax.axvline(x=sst_enso_n_p90, ymin=0, ymax=1, linewidth=1, color='b', linestyle='-.')
ax.axvline(x=sst_enso_n_mean, ymin=0, ymax=1, linewidth=1, color='b', linestyle='-.')
ax.legend()
ax.set_xlabel('SSTa', fontsize=14)
ax.set_ylabel('Probability density', fontsize=14)
ax.set_title(r'Histogram -- East PNG 5S 180E', fontsize=16, y=1.02) 
# of IQ: $\mu=100$, $\sigma=15$')
# Tweak spacing to prevent clipping of ylabel
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
fig.tight_layout()
plt.savefig(figfile_hist,bbox_inches='tight', format='eps', dpi=300)
#'''
#plt.show()




