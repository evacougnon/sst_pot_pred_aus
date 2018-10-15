import numpy as np
import pandas as pd
import xarray as xr

from datetime import date
from calendar import monthrange
import time as time

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_ClimateIndex as eac_CI


# specify if using the DPS outout or the satelite obs
WHAT = 'DPS'
#WHAT = 'OBS'

# load data
if WHAT == 'DPS':
    header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
    fname1 = header + 'NumDays_LaN_p90p10_Aus_dps.nc'
    fname2 = header + 'NumDays_ElN_p90p10_Aus_dps.nc'
    fname_sig_LaN = header + 'NumDays_LaN_p90p10_Aus_95sign_2500_dps'
    fname_sig_ElN = header + 'NumDays_ElN_p90p10_Aus_95sign_2500_dps'
    figfile = header + 'map_ssta_percentage_sig_dps_enso_Aus.png'
    lat_map = xr.open_dataset(fname1)['yt_ocean']
    lon_map = xr.open_dataset(fname1)['xt_ocean'] + 360 # to use with DPS model
elif WHAT == 'OBS':
    fname1='/home/ecougnon/ana/NumDays_LaN_p90p10_Aus_filter.nc'
    fname2='/home/ecougnon/ana/NumDays_ElN_p90p10_Aus_filter.nc'
    fname_sig_LaN = '/home/ecougnon/ana/NumDays_LaN_p90p10_Aus_95sign_2500_filter'
    fname_sig_ElN = '/home/ecougnon/ana/NumDays_ElN_p90p10_Aus_95sign_2500_filter'
    figfile = '/home/ecougnon/ana/map_ssta_percentage_enso_oz_sign_filter.png'
    lat_map = xr.open_dataset(fname1)['lat']
    lon_map = xr.open_dataset(fname1)['lon']

# percentage of days above/below p90/p10
p_p90_LaN = xr.open_dataset(fname1)['p_p90_LaN']
p_p10_LaN = xr.open_dataset(fname1)['p_p10_LaN']
p_p90_ElN = xr.open_dataset(fname2)['p_p90_ElN']
p_p10_ElN = xr.open_dataset(fname2)['p_p10_ElN']
# number of days above/below p90/p10
n_p90_LaN = xr.open_dataset(fname1)['n_p90_LaN']
n_p10_LaN = xr.open_dataset(fname1)['n_p10_LaN']
n_p90_ElN = xr.open_dataset(fname2)['n_p90_ElN']
n_p10_ElN = xr.open_dataset(fname2)['n_p10_ElN']
# limits for a 95% significance
data_LaN = np.load(fname_sig_LaN + '.npz')
data_ElN = np.load(fname_sig_ElN + '.npz')
rt_low_LaN_p90 = data_LaN['rt_low_LaN_p90']
rt_high_LaN_p90 = data_LaN['rt_high_LaN_p90']
rt_low_LaN_p10 = data_LaN['rt_low_LaN_p10']
rt_high_LaN_p10 = data_LaN['rt_high_LaN_p10']
rt_low_ElN_p90 = data_ElN['rt_low_ElN_p90']
rt_high_ElN_p90 = data_ElN['rt_high_ElN_p90']
rt_low_ElN_p10 = data_ElN['rt_low_ElN_p10']
rt_high_ElN_p10 = data_ElN['rt_high_ElN_p10']

mask_LaN_p90 = np.ones(rt_low_LaN_p90.shape)
mask_LaN_p90 = np.ma.masked_where((n_p90_LaN>rt_low_LaN_p90) & \
                                  (n_p90_LaN<rt_high_LaN_p90), \
                                  mask_LaN_p90)
mask_LaN_p10 = np.ones(rt_low_LaN_p10.shape)
mask_LaN_p10 = np.ma.masked_where((n_p10_LaN>rt_low_LaN_p10) & \
                                  (n_p10_LaN<rt_high_LaN_p10), \
                                  mask_LaN_p10)
mask_ElN_p90 = np.ones(rt_low_ElN_p90.shape)
mask_ElN_p90 = np.ma.masked_where((n_p90_ElN>rt_low_ElN_p90) & \
                                  (n_p90_ElN<rt_high_ElN_p90), \
                                  mask_ElN_p90)
mask_ElN_p10 = np.ones(rt_low_ElN_p10.shape)
mask_ElN_p10 = np.ma.masked_where((n_p10_ElN>rt_low_ElN_p10) & \
                                  (n_p10_ElN<rt_high_ElN_p10), \
                                  mask_ElN_p10)

if WHAT == 'OBS':
# usefull numbers
    MinYear = 2003 #1982
    MaxYear = 2017 #2016
    NumYears = MaxYear-MinYear+1
    MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                     # case the first year is a leap year
    NumDays = 365*NumYears + MaxNumLeapYear
# time vector for plotting
#    dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,12,31).toordinal()+1)
########################################
## CLIMATE MODE
########################################
######################################
# NINO34 index from the smoothed (1/4 degree to 1degree)
# NOAA OISSTV2 daily dataset
#########################################
    ds_oisst = xr.open_dataset('/home/ecougnon/ana/SST_smooth_1deg_NINO34region.nc')['__xarray_dataarray_variable__']. \
                  sel(lat=slice(-5,5), \
                      lon=slice(190,240), \
                      time=slice('2003-01-01','2017-10-31')). \
                  resample('1MS', dim='time', how='mean')
    clim_oisst = ds_oisst.groupby('time.month').mean('time')
    ssta_oisst = ds_oisst.groupby('time.month') - clim_oisst
    dssta_oisst = np.array(ssta_oisst)
#    dssta_oisst = np.apply_along_axis(signal.detrend,0,ssta_oisst, type='linear')
    nino34 = np.mean(dssta_oisst,axis=(1,2))
    mtim1e = pd.date_range('2003-01-01','2017-10-31',name='time',freq='M')
# warning finishes in Nov 2017 not Dec
    dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,10,31).toordinal()+1)
# make the index daily to be applied on the daily output
    nino34_d= np.nan*np.zeros(len(dtime))
    m=0
    y=0
    d=0
    for yy in np.arange(0,NumYears):
        if (yy==NumYears-1):
            for mm in np.arange(1,10+1):
                nino34_d[d:d+monthrange(MinYear+y,mm)[1]] = nino34[m]
                m = m + 1
                d = d + monthrange(MinYear+y,mm)[1]
        else:
            for mm in np.arange(1,12+1):
                nino34_d[d:d+monthrange(MinYear+y,mm)[1]] = nino34[m]
                m = m + 1
                d = d + monthrange(MinYear+y,mm)[1]
        y = y + 1
    enso_p_count = np.count_nonzero(nino34_d>=0.4)
    enso_n_count = np.count_nonzero(nino34_d<=-0.4)
    enso_neu_count = np.count_nonzero((nino34_d>-0.4) & (nino34_d<0.4))

    '''
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
    enso_p_count = np.count_nonzero(enso_daily>=0.75)
    enso_n_count = np.count_nonzero(enso_daily<=-0.75)
    enso_neu_count = np.count_nonzero((enso_daily>-0.75) & (enso_daily<0.75))
    '''

elif WHAT == 'DPS':
# usefull numbers
    MinYear = 2003
    MaxYear = 2017
    NumYear = MaxYear - MinYear+1
# warning finishes in Nov 2017 not Dec
    dtime = np.arange(date(MinYear,1,1).toordinal(),date(MaxYear,11,30). \
                      toordinal()+1)
    NumDays = len(dtime)
########################################
## CLIMATE MODE
########################################
# calc NINO3.4 index from DPS
###################################################
    nino34_id, nino34_ph = eac_CI.nino34_index_dps()
    mtime = pd.date_range('2003-01-01','2017-12-01',name='time',freq='M')
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
    enso_p_count = np.count_nonzero(nino34_d==1)
    enso_n_count = np.count_nonzero(nino34_d==-1)
    enso_neu_count = np.count_nonzero(nino34_d==0)

################################################
# plotting
################################################
domain = [-55, 90, 10, 180] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-55, 90, 10, 180]
dlat = 10
dlon = 30
if WHAT == 'OBS':
    llat, llon = np.meshgrid(lat_map, lon_map)
elif WHAT == 'DPS':
    llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'
bin_col = 0.1
bin_bar = 0.5

my_dpi = 300

ax=plt.figure(figsize=(4000/my_dpi,4000/my_dpi))
#ax=plt.figure(figsize=(15,17)) 
plt.clf()
plt.subplot(2,2,1, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, (n_p90_ElN/enso_p_count)*100, \
             levels=np.arange(0,60+1,10), cmap=plt.cm.Reds)
cb=plt.colorbar(ticks=np.arange(0,60+1,10),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.contourf(lonproj, latproj, mask_ElN_p90, hatches = '.', alpha=0)
plt.title('a) % of days above p90$_{ALL}$ during EN', \
          fontsize=16, y=1.02)

plt.subplot(2,2,3, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, (n_p10_ElN/enso_p_count)*100, \
             levels=np.arange(0,60+1,10), cmap=plt.cm.Blues)
cb=plt.colorbar(ticks=np.arange(0,60+1,10),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.contourf(lonproj, latproj, mask_ElN_p10, hatches = '.', alpha=0)
plt.title('b) % of days below p10$_{ALL}$ during EN', \
          fontsize=16, y=1.02)

plt.subplot(2,2,2, facecolor=bg_col)
#plt.figure()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, (n_p90_LaN/enso_n_count)*100, \
             levels=np.arange(0,60+1,10), cmap=plt.cm.Reds)
cb=plt.colorbar(ticks=np.arange(0,60+1,10),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.contourf(lonproj, latproj, mask_LaN_p90, hatches = '.', alpha=0)
plt.title('c) % of days above p90$_{ALL}$ during LN', \
          fontsize=16, y=1.02)

plt.subplot(2,2,4, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, (n_p10_LaN/enso_n_count)*100, \
             levels=np.arange(0,60+1,10), cmap=plt.cm.Blues)
cb=plt.colorbar(ticks=np.arange(0,60+1,10),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.contourf(lonproj, latproj, mask_LaN_p10, hatches = '.', alpha=0)
plt.title('d) % of days below p10$_{ALL}$ during LN', \
          fontsize=16, y=1.02)


plt.savefig(figfile,bbox_inches='tight', format='png', dpi=300)
plt.show()



