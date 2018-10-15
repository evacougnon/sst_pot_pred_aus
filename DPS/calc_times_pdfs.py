'''
    time series analyses for specific points chosen
    around Australia to compare to different climate
    mode -- starting with ENSO

	code based on ts_ana_pdfs.py used with the observations

    Author: Eva A Cougnon
    Created: Apr 2018
    Last modified:  

'''

# load required modules

import numpy as np
from scipy import io # read/load mat files
from scipy import signal  # detrend
from scipy.stats import pearsonr
import xarray as xr
import pandas as pd
from datetime import date
from calendar import monthrange

from matplotlib import pyplot as plt

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eric_oliver as eo
import eac_useful as eac
import eac_ClimateIndex as eac_CI


# calc NINO3.4 index
nino34_id, nino34_ph = eac_CI.nino34_index_dps()
mtime = pd.date_range('2003-01-01','2017-12-01',name='time',freq='M')

# define indexes for lat lon of one pixel
lat_px =  -30 #-5 #-30 #-42 #-30 # deg N
lon_px = 112 - 360 # 180 #112 #151 #112 # deg E

# fig files saving names
header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
figfile_ts = header + 'ts_ssta_enso_30-112.png'
#figfile_hist = header + 'hist_ssta_enso_42-151.png'

# usefull numbers
MinYear = 2003
MaxYear = 2017
NumYear = MaxYear - MinYear+1
# warning finishes in Nov 2017 not Dec
dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,11,30).toordinal()+1)

# make the index daily to be applied on the daily output
nino34_d= np.nan*np.zeros(len(dtime))
m=0
y=0
d=0
for yy in np.arange(0,NumYear):
    if (yy==NumYear-1):
        for mm in np.arange(1,11+1):
            nino34_d[d:d+monthrange(MinYear+y,mm)[1]] = nino34_ph[m]
            m = m + 1
            d = d + monthrange(MinYear+y,mm)[1]
    else:
        for mm in np.arange(1,12+1):
            nino34_d[d:d+monthrange(MinYear+y,mm)[1]] = nino34_ph[m]
            m = m + 1
            d = d + monthrange(MinYear+y,mm)[1]
    y = y + 1

# save indexes where enso is +ve (p) and -ve (n)
#nino34_p_id = np.nonzero(nino34_d>=0.4)
#nino34_n_id = np.nonzero(nino34_d<=-0.4)
#nino34_neu_id = np.nonzero((nino34_d>-0.4) & (nino34_d<0.4))


############################################################
# load the time series
############################################################
SSTa_d = xr.open_dataset(header + 'ssta_reana_ETKF_mem001_20032017_daily_Aus.nc') \
                        ['dsst_mdl'].sel(time=slice('2003-01-01','2017-11-30')). \
                                     sel(yt_ocean=lat_px, xt_ocean=lon_px, \
                                         method='nearest')
#SSTa_d = xr.open_dataset('/home/ecougnon/data/DPS/forecast/sst_daily_Aus_ensemble.nc')['sst']. \
#            sel(yt_ocean=lat_px, xt_ocean=lon_px, \
#                method='nearest'). \
#            sel(time=slice('2007-01-01','2010-12-31'))

'''
# monthly data
gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'
yt_ocean = xr.open_dataset(gname)['yt_ocean']
yt_ocean = yt_ocean.sel(yt_ocean=lat_px, method='nearest')
xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
xt_ocean = xt_ocean.sel(xt_ocean=lon_px, method='nearest')
SST_m = eac.read_netcdfs(header + 'ocean_daily_SST_enkf9_mem001_20??-20??.nc', \
                         dim='time', lat=yt_ocean, lon=xt_ocean, \
                         transform_func=lambda ds: ds.resample('1MS', \
                                                               dim='time', \
                                                               how='mean'))
SST_m = SST_m.squeeze('st_ocean')
SST_m = SST_m.sel(time=slice('2003-01-01','2017-11-01'))
# calc monthly climatology
clim_m = SST_m.groupby('time.month').mean('time')
ssta_m = SST_m.groupby('time.month') - clim_m
dssta_m = np.apply_along_axis(signal.detrend,0,ssta_m['temp'], type='linear')
'''

# calc stats from the point time serie
ssta_ts_mean = np.mean(SSTa_d)
ssta_ts_std = np.std(SSTa_d)
ssta_ts_p10 = np.percentile(SSTa_d,10)
ssta_ts_p90 = np.percentile(SSTa_d,90)
ssta_ts_med = np.median(SSTa_d)
ssta_ts_max = np.max(SSTa_d)
ssta_ts_min = np.min(SSTa_d)

# calc stats from the point time serie
# considering +ve or -ve ENSO
ssta_enso_n = SSTa_d.copy() # to use only on numpy array
ssta_enso_p = SSTa_d.copy()
ssta_enso_n[nino34_d==1]=np.nan
ssta_enso_n[nino34_d==0]=np.nan
ssta_enso_p[nino34_d==-1]=np.nan
ssta_enso_p[nino34_d==0]=np.nan

'''
ssta_enso_n[nino34_p_id[0][:]] = np.nan
ssta_enso_n[nino34_neu_id[0][:]] = np.nan
ssta_enso_p[nino34_n_id[0][:]] = np.nan
ssta_enso_p[nino34_neu_id[0][:]] = np.nan
'''
ssta_enso_p_mean = np.nanmean(ssta_enso_p)
ssta_enso_p_std = np.nanstd(ssta_enso_p)
ssta_enso_p_p10 = np.nanpercentile(ssta_enso_p,10)
ssta_enso_p_p90 = np.nanpercentile(ssta_enso_p,90)
ssta_enso_p_med = np.nanmedian(ssta_enso_p)
ssta_enso_p_max = np.nanmax(ssta_enso_p)
ssta_enso_p_min = np.nanmin(ssta_enso_p)

ssta_enso_n_mean = np.nanmean(ssta_enso_n)
ssta_enso_n_std = np.nanstd(ssta_enso_n)
ssta_enso_n_p10 = np.nanpercentile(ssta_enso_n,10)
ssta_enso_n_p90 = np.nanpercentile(ssta_enso_n,90)
ssta_enso_n_med = np.nanmedian(ssta_enso_n)
ssta_enso_n_max = np.nanmax(ssta_enso_n)
ssta_enso_n_min = np.nanmin(ssta_enso_n)

#'''
#############################################
## PLOTTING
############################################
dtime_ = pd.date_range('2003-01-01','2017-11-30',name='time',freq='D')

my_dpi = 300

plt.figure(figsize=(4000/my_dpi,3000/my_dpi))
plt.clf()
# plot p10 threshold
ax1=plt.subplot(2,1,1)
plt.plot(mtime, nino34_id, color='k')
ax1.set_xlim([mtime[0], mtime[-1]])
ax1.set_ylim([-3, 4.5])
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
ax1.fill_between(mtime, nino34_id, where=nino34_ph==1, facecolor='red', alpha=0.5, \
                 interpolate=True)
ax1.fill_between(mtime, nino34_id, where=nino34_ph==-1, facecolor='blue', alpha=0.5, \
                 interpolate=True)
plt.title('NINO3.4 Index', fontsize=16, y=1.02)
plt.grid()
#plt.gcf().autofmt_xdate() # make the x axis more readable 
ax2=plt.subplot(2,1,2)
plt.plot(dtime_,SSTa_d, color='0.7')
plt.plot(dtime_,ssta_enso_n, color='b')
plt.plot(dtime_,ssta_enso_p, color='r')
plt.axhline(y=ssta_ts_p10, color='k', linestyle=':', linewidth=1, \
            label='time series')
plt.axhline(y=ssta_enso_n_p10, color='b', linestyle=':', linewidth=1, \
            label='-ENSO')
plt.axhline(y=ssta_enso_p_p10, color='r', linestyle=':', linewidth=1, \
            label='+ENSO')
plt.axhline(y=ssta_ts_p90, color='k', linestyle=':', linewidth=1)
plt.axhline(y=ssta_enso_n_p90, color='b', linestyle=':', linewidth=1)
plt.axhline(y=ssta_enso_p_p90, color='r', linestyle=':', linewidth=1)
#ax2.legend()
plt.title('SSTa at 30S 112E -- WA', fontsize=16, y=1.02)
ax2.set_xlim([dtime_[0], dtime_[-1]])
ax2.set_ylim([-4, 4])
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
#plt.grid()
plt.savefig(figfile_ts,bbox_inches='tight', format='png', dpi=300)
#'''

#################
# histogram
#################
num_bins = 100 
fig, ax = plt.subplots()
# the histogram of the data
#n, bins, patches = ax.hist([sst_enso_n[~np.isnan(sst_enso_n)]], num_bins, \
#normed=1, histtype='bar',stacked=True)
n,bins,patches = plt.hist(ssta_enso_n[~np.isnan(ssta_enso_n)], num_bins, \
                          normed = True, histtype='step', linewidth=2, \
                          color='b', alpha=0.7, stacked=True, label='-ve ENSO')
n,bins,patches = plt.hist(ssta_enso_p[~np.isnan(ssta_enso_p)], num_bins, \
                          normed = True, histtype='step', linewidth=2, \
                          color='r', alpha=0.7, stacked=True, label='+ve ENSO')
n,bins,patches = plt.hist(SSTa_d, num_bins, normed = True, histtype='step', \
                          linewidth=2, color='k', alpha=0.7, \
                          stacked=True, label='timeseries')
# add a 'best fit' line
#y = mlab.normpdf(bins, sst_ts_mean, sst_ts_std)
#ax.plot(bins, y, '--')
ax.legend()
ax.axvline(x=ssta_ts_p10, ymin=0, ymax=1, linewidth=1, color='k', linestyle=':')
ax.axvline(x=ssta_ts_p90, ymin=0, ymax=1, linewidth=1, color='k', linestyle=':')
ax.axvline(x=ssta_ts_mean, ymin=0, ymax=1, linewidth=1, color='k', linestyle=':')
ax.axvline(x=ssta_enso_p_p10, ymin=0, ymax=1, linewidth=1, color='r', \
           linestyle=':')
ax.axvline(x=ssta_enso_p_p90, ymin=0, ymax=1, linewidth=1, color='r', \
           linestyle=':')
ax.axvline(x=ssta_enso_p_mean, ymin=0, ymax=1, linewidth=1, color='r', \
           linestyle=':')
ax.axvline(x=ssta_enso_n_p10, ymin=0, ymax=1, linewidth=1, color='b', \
           linestyle=':')
ax.axvline(x=ssta_enso_n_p90, ymin=0, ymax=1, linewidth=1, color='b', \
           linestyle=':')
ax.axvline(x=ssta_enso_n_mean, ymin=0, ymax=1, linewidth=1, color='b', \
           linestyle=':')
ax.set_xlabel('SSTa', fontsize=14)
ax.set_ylabel('Probability density', fontsize=14)
ax.set_title(r'Histogram -- TAS 42S 151E', fontsize=16, y=1.02)
# of IQ: $\mu=100$, $\sigma=15$')
# Tweak spacing to prevent clipping of ylabel
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
fig.tight_layout()
#plt.savefig(figfile_hist,bbox_inches='tight', format='png', dpi=300)
plt.show()



