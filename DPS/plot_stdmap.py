'''

    plot maps of STD from the monthly SSTa and the monthly
	p10 and p90 based on the daily SSTa

'''

# load required modules

import numpy as np
import pandas as pd
import xarray as xr

from datetime import date
from calendar import monthrange
import time as time

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
figfile = header +'SSTa_monthly_std_stats_dps.eps'
###############################
# Load data
################################
# define the region
lat_min = -55
lat_max = 10
# !!! LONGITUDE in the model!!! from -280 to 80 deg E
lon_min = 90 - 360
lon_max = 180 - 360
# grid info
gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'
yt_ocean = xr.open_dataset(gname)['yt_ocean']
yt_ocean = yt_ocean.sel(yt_ocean=slice(lat_min,lat_max))
xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
xt_ocean = xt_ocean.sel(xt_ocean=slice(lon_min,lon_max))
# usefull numbers
MinYear = 2003
MaxYear = 2017
NumYear = MaxYear - MinYear+1
# warning finishes in Nov 2017 not Dec
dtime = pd.date_range('2003-01-01','2017-11-30',name='time',freq='D')
NumDays = len(dtime)
month_vec = pd.date_range('2003-01-01','2017-11-30',name='time',freq='M')
str_id=np.empty(len(month_vec))
end_id=np.empty(len(month_vec))
k=0
for yy in range(MinYear,MaxYear+1):
    if yy==MaxYear:
        for mm in range(1,11+1):
            tmp1 = np.where(dtime == np.datetime64(date(yy,mm,1)))
            str_id[k] = tmp1[0][0]
            tmp2 = np.where(dtime == np.datetime64(date(yy,mm, \
                                                        monthrange(yy,mm)[1])))
            end_id[k] = tmp2[0][0]
            k = k +1
    else:
        for mm in range(1,12+1):
            tmp1 = np.where(dtime == np.datetime64(date(yy,mm,1)))
            str_id[k] = tmp1[0][0]
            tmp2 = np.where(dtime == np.datetime64(date(yy,mm, \
                                                        monthrange(yy,mm)[1])))
            end_id[k] = tmp2[0][0]
            k = k +1
###########################
# Daily SSTa
###########################
SSTa_d = xr.open_dataset(header + 'ssta_reana_ETKF_mem001_20032017_daily_Aus.nc') \
                        ['dsst_mdl'].sel(time=slice('2003-01-01','2017-11-30'), \
                                         yt_ocean=slice(lat_min,lat_max), \
                                         xt_ocean=slice(lon_min,lon_max))

# calc the monthly stats (mean, p10, p90)
SSTa_m_mean = SSTa_d.resample('1MS', dim='time', how='mean')
'''
# only output the 12 month value. Does the percentile on all Jans, Febs, ...
SSTa_m_p10 = SSTa_d.groupby('time.month').reduce(np.nanpercentile, dim='time', q=10)
'''
SSTa_m_p10 = np.nan*np.zeros((len(month_vec),len(yt_ocean), len(xt_ocean)))
SSTa_m_p90 = np.nan*np.zeros((len(month_vec),len(yt_ocean), len(xt_ocean)))
for mm in range(0,len(month_vec)):
    SSTa_m_p10[mm,:,:] = SSTa_d[int(str_id[mm]):int(end_id[mm]+1),:,:] \
                         .quantile(0.1,dim='time')
    SSTa_m_p90[mm,:,:] = SSTa_d[int(str_id[mm]):int(end_id[mm]+1),:,:] \
                         .quantile(0.9,dim='time')

####################################
# plotting
####################################
# setting
domain = [-55, 90, 10, 180] #[-80, -180, 85, 180] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(xt_ocean + 360,yt_ocean)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure(figsize=(15,5)) #(19,7))
plt.clf()

plt.subplot(1,2,1, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, np.nanstd(SSTa_d,axis=0), \
             levels=np.arange(0,1.2+0.1,0.1), cmap=plt.cm.YlOrBr)
cb=plt.colorbar(ticks=np.arange(0,1.2+0.1,0.2),shrink=0.8)
cb.ax.tick_params(labelsize=14)
plt.title('STD monthly SSTa -- DPS enkf-9', \
          fontsize=14, y=1.02)

plt.subplot(1,2,2, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, \
             np.nanstd(SSTa_m_p90,axis=0) - np.nanstd(SSTa_m_p10,axis=0), \
             levels=np.arange(-0.2,0.2+0.05,0.05), \
#(0,1.2+0.1,0.1), \
             cmap=plt.cm.bwr) #YlOrBr)
#cb=plt.colorbar(ticks=np.arange(0,1.2+0.1,0.2),shrink=0.5)
cb=plt.colorbar(ticks=np.arange(-0.2,0.2+0.1,0.1),shrink=0.8)
cb.ax.tick_params(labelsize=14)
plt.title('STD monthly p90 SSTa - STD monthly p10 SSTa', \
          fontsize=14, y=1.02)

plt.savefig(figfile,bbox_inches='tight', format='eps', dpi=300)
plt.show()



