'''
    read DPS monthly model output from the reanalysis run
    and compar eto the NOAA OISST V2 used in the potential 
    predictability analysis

    ONLY for the Australian region (-55:10N -- 90:180E)

   Author: Eva A Cougnon
   Created: Oct 2017
   Last modification:  

'''

# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import datetime
from eofs.xarray import Eof
from scipy import signal # detrend
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm


# load data
fname_mdl = '/home/ecougnon/data/DPS/reanalysis/ssta_reana_20022016_daily_Aus.nc'
lat = xr.open_dataset(fname_mdl)['yt_ocean']
lon = xr.open_dataset(fname_mdl)['xt_ocean']+360 
tim = xr.open_dataset(fname_mdl)['time']
sst_mdl = xr.open_dataset(fname_mdl)['dsst_mdl']
'''
for la in range(0,len(lat)):
    for lo in range(0,len(lon)):
        valid = ~np.isnan(sst_mdl[:,la,lo])
        if (valid.any()==True):
            sst_mdl[:,la,lo] = signal.detrend(sst_mdl[:,la,lo], axis=0, \
                                              type='linear')
        elif (valid.all()==False):
            sst_mdl[:,la,lo] = np.nan
'''

fname_obs = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
#SST_monthly_extremes_Aus'
# the file is nearly 2gb!!!.... takes a while to load!
data = np.load(fname_obs + '.npz')
lon_obs = data['lon_map']
lat_obs = data['lat_map']
sst_obs_X = data['SST'].item()
# use only the monthly mean
tim_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
tim_str = np.where(tim_vec == np.datetime64(datetime.datetime(2002,1,31)))
tim_end = np.where(tim_vec == np.datetime64(datetime.datetime(2016,6,30)))
sst_obs_tmp = sst_obs_X['TMM'][:,:,tim_str[0][0]:tim_end[0][0]+1]
# swap axes to get time, lat, lon
sst_obs_tmp = np.swapaxes(np.swapaxes(sst_obs_tmp,2,0),1,2)
# change to a xarray to apply the same
X = len(lon_obs)
Y = len(lat_obs)
nm = len(tim)
SST_OBS = xr.Dataset({'sst_obs':(('time','lat','lon'), np.zeros((nm,Y,X)))}, \
                     {'time': tim, 'lat': lat_obs, 'lon': lon_obs})
SST_OBS['sst_obs'][:,:,:]=sst_obs_tmp[:,:,:]
sst_obs=xr.DataArray(SST_OBS['sst_obs'])
'''
for la in range(0,len(lat_obs)):
    for lo in range(0,len(lon_obs)):
        valid = ~np.isnan(sst_obs[:,la,lo])
        if (valid.any()==True):
            sst_obs[:,la,lo] = signal.detrend(sst_obs[:,la,lo], axis=0, \
                                               type='linear')
        elif (valid.all()==False):
            sst_obs[:,la,lo] = np.nan
'''

# EOF for model
coslat_mdl = np.cos(np.deg2rad(sst_mdl.coords['lat'].values))
wgts_mdl = np.sqrt(coslat_mdl)[..., np.newaxis]
solver_mdl = Eof(sst_mdl, weights=wgts_mdl, center=True)
lambdas_mdl=solver_mdl.eigenvalues()
vf_mdl = solver_mdl.varianceFraction()
Nerror_mdl = solver_mdl.northTest(vfscaled=True)
pcs_mdl = solver_mdl.pcs() #(time, mode)
eofs_mdl = solver_mdl.eofs()
# EOF for obs
coslat_obs = np.cos(np.deg2rad(sst_obs.coords['lat'].values))
wgts_obs = np.sqrt(coslat_obs)[..., np.newaxis]
solver_obs = Eof(sst_obs, weights=wgts_obs, center=True)
lambdas_obs=solver_obs.eigenvalues()
vf_obs = solver_obs.varianceFraction()
Nerror_obs = solver_obs.northTest(vfscaled=True)
pcs_obs = solver_obs.pcs() #(time, mode)
eofs_obs = solver_obs.eofs()




## plotting
domain = [-55, 90, 10, 180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
domain_draw = [-55, 90, 10, 180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon_obs, llat_obs = np.meshgrid(lon_obs, lat_obs)
llon_mdl, llat_mdl = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'
lev = np.hstack((np.arange(-0.06,-0.0005+0.0005,0.001), \
                 np.arange(0.0005,0.06+0.0005,0.001)))

# PC1
plt.figure
pcs_mdl[:, 0].plot(color='b', linewidth=2)
pcs_obs[:, 0].plot(color='k', linewidth=2)
ax = plt.gca()
ax.axhline(0, color='k')
#ax.set_ylim(-3, 3)
ax.set_xlabel('Year')
ax.legend(['obs','mdl'])
ax.set_ylabel('PC amplitude00')
ax.set_title('PC1 Time Series', fontsize=16)
# Eigenvalue spectrum
plt.figure()
plt.errorbar(np.arange(0,10),vf_obs[:10]*100, \
             yerr=Nerror_obs[:10]*100,fmt='.')
plt.errorbar(np.arange(0,10),vf_mdl[:10]*100, \
             yerr=Nerror_mdl[:10]*100,fmt='.')
plt.title('Eigenvalues spectrum -- North et al. 1982')
plt.legend(['obs','mdl'])
plt.xlabel('EOF modes')
plt.ylabel('Eigenvalue (%)')
plt.grid()
#plt.show()
# EOF map
plt.figure(figsize=(11,11)) #12,11)) #7,11)) #(12,6)) # (10,8)
plt.clf()
ax1 = plt.subplot(2,2,1)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_mdl, llat_mdl)
plt.contourf(lonproj, latproj, eofs_mdl[2,:,:], \
             levels= lev, \
             cmap=plt.cm.seismic)
#cd = plt.colorbar()
cb=plt.colorbar(ticks=lev,shrink=0.9)
#cb.ax1.tick_params(labelsize=14)
plt.title('EOF 3, mdl', \
          fontsize=16, y=1.04)
ax2 = plt.subplot(2,2,2)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_mdl, llat_mdl)
plt.contourf(lonproj, latproj, eofs_mdl[3,:,:], levels=lev,  cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev, shrink=0.9)
#cb.ax2.tick_params(labelsize=14)
plt.title('EOF 4, mdl', \
          fontsize=16, y=1.04)
ax3 = plt.subplot(2,2,3)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_obs, llat_obs)
plt.contourf(lonproj, latproj, eofs_obs[2,:,:], \
             levels= lev, \
             cmap=plt.cm.seismic)
#cd = plt.colorbar()
cb=plt.colorbar(ticks=lev,shrink=0.9)
#cb.ax1.tick_params(labelsize=14)
plt.title('EOF 3, obs', \
          fontsize=16, y=1.04)
ax4 = plt.subplot(2,2,4)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_obs, llat_obs)
plt.contourf(lonproj, latproj, eofs_obs[3,:,:], levels=lev,  cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev, shrink=0.9)
#cb.ax2.tick_params(labelsize=14)
plt.title('EOF 4, obs', \
          fontsize=16, y=1.04)

#plt.savefig(figfile, bbox_inches='tight',dpi=300)
plt.show()



