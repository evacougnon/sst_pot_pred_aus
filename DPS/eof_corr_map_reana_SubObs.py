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
from eofs.xarray import Eof
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays
from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

figfile = '/home/ecougnon/data/DPS/reanalysis/CorrMap_mdl_obs_ssta_Aus.png'

# Load data
fname = '/home/ecougnon/data/DPS/reanalysis/EOFsolvers_mdl_obs_ssta_Aus'
data = np.load(fname + '.npz')
lat = data['lat_mdl']
lon = data['lon_mdl']
tim = data['tim_mdl']
pcs_mdl = data['pcs_mdl']
eofs_mdl = data['eofs_mdl']
pcs_obs = data['pcs_obs']
eofs_obs = data['eofs_obs']

t = len(tim)
X = len(lon)
Y = len(lat)

mode1_mdl=np.empty((t,Y,X))
mode2_mdl=np.empty((t,Y,X))
mode3_mdl=np.empty((t,Y,X))
mode4_mdl=np.empty((t,Y,X))
mode1_mdl.fill(np.nan)
mode2_mdl.fill(np.nan)
mode3_mdl.fill(np.nan)
mode1_mdl.fill(np.nan)
mode1_obs=np.empty((t,Y,X))
mode2_obs=np.empty((t,Y,X))
mode3_obs=np.empty((t,Y,X))
mode4_obs=np.empty((t,Y,X))
mode1_obs.fill(np.nan)
mode2_obs.fill(np.nan)
mode3_obs.fill(np.nan)
mode1_obs.fill(np.nan)
## spatial correlation maps of the 4 leading modes
for i in range(0,t):
    mode1_mdl[i,:,:] = pcs_mdl[i, 0] * eofs_mdl[0,:,:]
    mode2_mdl[i,:,:] = pcs_mdl[i, 1] * eofs_mdl[1,:,:]
    mode3_mdl[i,:,:] = pcs_mdl[i, 2] * eofs_mdl[2,:,:]
    mode4_mdl[i,:,:] = pcs_mdl[i, 3] * eofs_mdl[3,:,:]
    mode1_obs[i,:,:] = pcs_obs[i, 0] * eofs_obs[0,:,:]
    mode2_obs[i,:,:] = pcs_obs[i, 1] * eofs_obs[1,:,:]
    mode3_obs[i,:,:] = pcs_obs[i, 2] * eofs_obs[2,:,:]
    mode4_obs[i,:,:] = pcs_obs[i, 3] * eofs_obs[3,:,:]
# all location need to be together in order to apply st.spearmanr function
# I tried to used np.apply_along_axis but it wouldn' take the correct axis 
# for the second array argument
mode1_mdl = mode1_mdl.reshape((t,Y*X))
mode2_mdl = mode2_mdl.reshape((t,Y*X))
#stack(z=('yt_ocean','xt_ocean'))
mode3_mdl = mode3_mdl.reshape((t,Y*X))
mode4_mdl = mode4_mdl.reshape((t,Y*X))
mode1_obs = mode1_obs.reshape((t,Y*X))
mode2_obs = mode2_obs.reshape((t,Y*X))
mode3_obs = mode3_obs.reshape((t,Y*X))
mode4_obs = mode4_obs.reshape((t,Y*X))

corr_map_mode1 = np.empty(X*Y)
corr_map_mode2 = np.empty(X*Y)
corr_map_mode3 = np.empty(X*Y)
corr_map_mode4 = np.empty(X*Y)
corr_map_mode1.fill(np.nan)
corr_map_mode2.fill(np.nan)
corr_map_mode3.fill(np.nan)
corr_map_mode4.fill(np.nan)
for ii in range(0,(X*Y)):
    corr_map_mode1[ii], tmp = st.pearsonr(mode1_mdl[:,ii], mode1_obs[:,ii]) 
    corr_map_mode2[ii], tmp = st.pearsonr(mode2_mdl[:,ii], mode2_obs[:,ii])
    corr_map_mode3[ii], tmp = st.pearsonr(mode3_mdl[:,ii], mode3_obs[:,ii])
    corr_map_mode4[ii], tmp = st.pearsonr(mode4_mdl[:,ii], mode4_obs[:,ii])
# change shape back to lat/lon
corr_map_mode1 = np.reshape(corr_map_mode1,(Y,X))
corr_map_mode2 = np.reshape(corr_map_mode2,(Y,X))
corr_map_mode3 = np.reshape(corr_map_mode3,(Y,X))
corr_map_mode4 = np.reshape(corr_map_mode4,(Y,X))

## plotting
domain = [-55, 90, 10, 180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
domain_draw = [-55, 90, 10, 180] #[-55, -270, 10, -180] #[-55, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon_obs, llat_obs = np.meshgrid(lon, lat)
llon_mdl, llat_mdl = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'
lev = np.hstack((np.arange(-1,-0.1+0.05,0.1), \
                 np.arange(0.1,1+0.05,0.1)))

plt.figure(figsize=(11,11)) 
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
plt.contourf(lonproj, latproj, corr_map_mode1, \
             levels= lev, \
             cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.9)
plt.title('Pearson Correlation map for EOF1', \
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
plt.contourf(lonproj, latproj, corr_map_mode2, levels=lev,  cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev, shrink=0.9)
plt.title('Pearson Correlation map for EOF 2', \
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
plt.contourf(lonproj, latproj, corr_map_mode3, levels= lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.9)
plt.title('Pearson Correlation map for EOF 3', \
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
plt.contourf(lonproj, latproj, corr_map_mode4, levels=lev,  cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev, shrink=0.9)
plt.title('Pearson Correlation map for EOF 4', \
          fontsize=16, y=1.04)

#plt.savefig(figfile, bbox_inches='tight',dpi=300)
plt.show()










## Case studies


