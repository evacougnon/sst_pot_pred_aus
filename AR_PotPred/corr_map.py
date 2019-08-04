'''
    plot the correlation maps for each EOF pattern/PCs with a 
    climate driver time series
'''

# load required modules

import numpy as np
import xarray as xr
import pandas as pd
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays
from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_ClimateIndex as eac_CI

# load data
header = '/home/ecougnon/ana/PotPred/EOF/'
fname = header + 'cov_decomp_1yr_ZFL_EOF_Aus_1deg_TMM'
data = np.load(fname + '.npz')
lat = data['lat']
lon = data['lon']
tim_vec = data['time_ts']
pcs = data['PC'].item()
eofs = data['MODES_MAP'].item()

var_PC = pcs['PC1_pred'].copy() 
var_EOF = eofs['pred_modes_map'].copy()
ind = 0 # index of the mode: mode 1 is 0

# useful dim
T = len(var_PC)
X = len(var_EOF[0,:,0])
Y = len(var_EOF[0,0,:])

## CLIMATE MODE
mei_monthly, mei_monthly_std = eac_CI.calc_month_index('MEI',1982,2016)
nino34_monthly, nino34_monthly_std = eac_CI.calc_month_index('NINO34',1982,2016)
dmi_monthly, dmi_monthly_std = eac_CI.calc_month_index('DMI',1982,2016)
sam_monthly, sam_monthly_std = eac_CI.calc_month_index('SAM',1982,2016)
nni_monthly, nni_monthly_std = eac_CI.calc_month_index('NNI',1982,2016)
eac_monthly, eac_str, eac_end = eac_CI.calc_month_eac(tim_vec)

# PC*EOF
mode_map = np.empty((T, X, Y))
mode_map.fill(np.nan)
corr_map_mei = np.empty(X*Y)
corr_map_mei.fill(np.nan)
corr_map_nino34 = np.empty(X*Y)
corr_map_nino34.fill(np.nan)
corr_map_dmi = np.empty(X*Y)
corr_map_dmi.fill(np.nan)
corr_map_sam = np.empty(X*Y)
corr_map_sam.fill(np.nan)
corr_map_nni = np.empty(X*Y)
corr_map_nni.fill(np.nan)
corr_map_eac = np.empty(X*Y)
corr_map_eac.fill(np.nan)
for tt in np.arange(0,T):
    mode_map[tt,:,:] = np.squeeze(var_EOF[ind,:,:]) * var_PC[tt]
mode_map_ = mode_map.reshape(T,X*Y)

for ll in np.arange(0,X*Y):
    corr_map_mei[ll], tmp = st.spearmanr(mode_map_[:,ll],mei_monthly)
    corr_map_nino34[ll], tmp = st.spearmanr(mode_map_[:,ll],nino34_monthly)
    corr_map_dmi[ll], tmp = st.spearmanr(mode_map_[:,ll],dmi_monthly)
    corr_map_sam[ll], tmp = st.spearmanr(mode_map_[:,ll],sam_monthly)
    corr_map_nni[ll], tmp = st.spearmanr(mode_map_[:,ll],nni_monthly)
#    corr_map_eac[ll], tmp = st.spearmanr(mode_map_[:,ll],eac_monthly)

corr_map_mei = corr_map_mei.reshape(X,Y)
corr_map_nino34 = corr_map_nino34.reshape(X,Y)
corr_map_dmi = corr_map_dmi.reshape(X,Y)
corr_map_sam = corr_map_sam.reshape(X,Y)
corr_map_nni = corr_map_nni.reshape(X,Y)
#corr_map_eac = corr_map_eac.reshape(X,Y)

# plot setting
domain = [-55, 90, 10, 180] #[-80, 0, 85, 360] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'

var_plot = corr_map_mei.copy()

ax=plt.figure(figsize=(11,11))
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, var_plot, levels=np.arange(-1,1+0.1,0.1), \
             cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=np.arange(-1,1+0.1,0.1),shrink=0.9)
cb.ax.tick_params(labelsize=14)
plt.title('correlation between mode 1 and MEI', \
          fontsize=16, y=1.08)
#plt.savefig(figfile,bbox_inches='tight', format='eps', dpi=300)
plt.show()




