'''
	calc EOF on the SSTa from the reanalysis and compares
	it with the AVHRR OISST V2 observations -- subset to match
	the model resolution (1deg in longitude, keeps changing in
	latitude)

   Author: Eva A Cougnon
   Created: Oct 2017
   Last modification:  

'''
# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import time as time
from eofs.xarray import Eof
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays
from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

figfile = '/home/ecougnon/data/DPS/reanalysis/ETKF/EOF3-4_ssta_Aus_20032016.png'

'''
# EOF solvers take about 5 minutes for the mdl and 26 minutes for the obs!!!
outfile = '/home/ecougnon/data/DPS/reanalysis/ETKF/EOFsolvers_mdl_obs_ssta_ETKF_Aus_20032016'

# Load data
fname_mdl = '/home/ecougnon/data/DPS/reanalysis/ETKF/ssta_reana_ETKF_mem001_20032017_daily_Aus.nc'
lat_mdl = xr.open_dataset(fname_mdl)['yt_ocean']
lon_mdl = xr.open_dataset(fname_mdl)['xt_ocean']+360
tim_mdl = xr.open_dataset(fname_mdl)['time']
tim_mdl = tim_mdl.sel(time=slice('2003-01-01', '2016-12-31'))
sst_mdl = xr.open_dataset(fname_mdl)['dsst_mdl']
sst_mdl = sst_mdl.sel(time=tim_mdl)

fname_obs = '/home/ecougnon/ana/SSTa_daily_Aus_20032016Dec.nc' #shift_time_dim.nc'
lat_obs = xr.open_dataset(fname_obs)['lat']
lat_obs = lat_obs.sel(lat=lat_mdl, method='nearest')
lon_obs = xr.open_dataset(fname_obs)['lon']
lon_obs = lon_obs.sel(lon=lon_mdl, method='nearest')
tim_obs = xr.open_dataset(fname_obs)['time']
#tim_obs = tim_obs.sel(time=slice('2002-01-01','2016-06-30'))
sst_obs = xr.open_dataset(fname_obs)['SSTa']
sst_obs = sst_obs.sel(time=tim_obs, lat=lat_obs,lon=lon_obs)
##taking ages!!!
# fname_obs = '/home/ecougnon/ana/SSTa_daily_Aus.nc'
#...
sst_obs = sst_obs.transpose('time','lat','lon')
sst_obs.to_netcdf('/home/ecougnon/ana/SSTa_daily_Aus_20032016Dec_shift_time_dim.nc')
################
t = len(tim_obs)
X = len(lon_obs)
Y = len(lat_obs)

## EOF
# EOF for mdl
t_key = time.time()
coslat_mdl = np.cos(np.deg2rad(sst_mdl.coords['yt_ocean'].values))
wgts_mdl = np.sqrt(coslat_mdl)[..., np.newaxis]
solver_mdl = Eof(sst_mdl, weights=wgts_mdl, center=True)
lambdas_mdl=solver_mdl.eigenvalues()
vf_mdl = solver_mdl.varianceFraction()
Nerror_mdl = solver_mdl.northTest(vfscaled=True)
pcs_mdl = solver_mdl.pcs() #(time, mode)
eofs_mdl = solver_mdl.eofs()
elapsed_key = time.time() - t_key
print('elapsed time for mdl EOF solver:', elapsed_key)
# EOF for obs
t_key = time.time()
coslat_obs = np.cos(np.deg2rad(sst_obs.coords['lat'].values))
wgts_obs = np.sqrt(coslat_obs)[..., np.newaxis]
solver_obs = Eof(sst_obs, weights=wgts_obs, center=True)
lambdas_obs=solver_obs.eigenvalues()
vf_obs = solver_obs.varianceFraction()
Nerror_obs = solver_obs.northTest(vfscaled=True)
pcs_obs = solver_obs.pcs() #(time, mode)
eofs_obs = solver_obs.eofs()
elapsed_key = time.time() - t_key
print('elapsed time for obs EOF solver:', elapsed_key)

## saving
np.savez(outfile, lambdas_mdl=lambdas_mdl, vf_mdl=vf_mdl, Nerror_mdl=Nerror_mdl, \
         pcs_mdl=pcs_mdl, eofs_mdl=eofs_mdl, lambdas_obs=lambdas_obs, \
         vf_obs=vf_obs, Nerror_obs=Nerror_obs, pcs_obs=pcs_obs, eofs_obs=eofs_obs, \
         lat_mdl=lat_mdl, lon_mdl=lon_mdl, tim_mdl=tim_mdl, lat_obs=lat_obs, \
         lon_obs=lon_obs, tim_obs=tim_obs)
'''

# Load data
fname = '/home/ecougnon/data/DPS/reanalysis/ETKF/EOFsolvers_mdl_obs_ssta_ETKF_Aus_20032016'
data = np.load(fname + '.npz')
lat = data['lat_mdl']
lon = data['lon_mdl']
tim = data['tim_mdl']
pcs_mdl = data['pcs_mdl']
eofs_mdl = data['eofs_mdl']
vf_mdl = data['vf_mdl']
Nerror_mdl = data['Nerror_mdl']
pcs_obs = data['pcs_obs']
eofs_obs = data['eofs_obs']
vf_obs = data['vf_obs']
Nerror_obs = data['Nerror_obs']



t = len(tim)
X = len(lon)
Y = len(lat)
'''

'''
# plotting
domain = [-55, 90, 10, 180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
domain_draw = [-55, 90, 10, 180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon_obs, llat_obs = np.meshgrid(lon, lat)
llon_mdl, llat_mdl = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'
lev = np.hstack((np.arange(-0.06,-0.0005+0.0005,0.001), \
                 np.arange(0.0005,0.06+0.0005,0.001)))

# PC1
pc = 3
SpearC, tmp = st.spearmanr(pcs_mdl[:, pc], pcs_obs[:, pc])
plt.figure
plt.plot(pcs_mdl[:, pc], color='b', linewidth=2)
plt.plot(pcs_obs[:, pc], color='k', linewidth=2)
ax = plt.gca()
ax.axhline(0, color='k')
#ax.set_ylim(-3, 3)
ax.set_xlabel('Year')
ax.legend(['mdl','obs'])
ax.set_ylabel('PC amplitude00')
ax.set_title('PC4 Time Series', fontsize=16)
plt.text(0.3, 0.1, 'Spearman Correlation coefficient:' + str(round(SpearC,3)), \
         ha='center', va='center', transform=ax.transAxes, \
         fontsize=14)
plt.grid()
plt.show()

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
plt.contourf(lonproj, latproj, eofs_mdl[0,:,:], \
             levels= lev, \
             cmap=plt.cm.seismic)
cb=plt.colorbar() #(ticks=lev,shrink=0.9)
plt.title('EOF 1, mdl', \
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
plt.contourf(lonproj, latproj, eofs_mdl[1,:,:], levels=lev,  cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev, shrink=0.9)
plt.title('EOF 2, mdl', \
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
plt.contourf(lonproj, latproj, eofs_obs[0,:,:], \
             levels= lev, \
             cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.9)
plt.title('EOF 1, obs', \
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
plt.contourf(lonproj, latproj, eofs_obs[1,:,:], levels=lev,  cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev, shrink=0.9)
plt.title('EOF 2, obs', \
          fontsize=16, y=1.04)

plt.savefig(figfile, bbox_inches='tight',dpi=300)
plt.show()






