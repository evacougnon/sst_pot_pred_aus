'''

    read and plot outputs from the potential predictability 
    analysis --> extract the diagonal of the covariance matrixes

'''


# import libraries
import numpy as np
import xarray as xr

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

from matplotlib.mlab import griddata
 
def grid(x, y, z, resX=100, resY=100):
    "Convert 3 column data to matplotlib grid"
    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z


fname = '../../ana/PotPred/PotPred_CovMtx_Aus_daily_1deg'
data = np.load(fname + '.npz')
lat = data['lat_NoLand']
lon = data['lon_NoLand']
cov_mtx_chunk = data['cov_mtx_chunk']
cov_noise = data['cov_noise']

#figfile = '/home/ecougnon/ana/PotPred/PPratio_SSTa_TMM_1yr_ZFL_sign_vSZ.png'
#PPratio_OTE_p90_1yr_ZFL_Aus.png'

######################################
# extract the diagonal for ppr calc
######################################
diag_tot = np.diag(cov_mtx_chunk)
diag_unp = np.diag(cov_noise)
diag_pred = np.diag(cov_mtx_chunk - cov_noise)

#######################################
# var to plot
######################################
var = diag_pred/diag_tot

##########################################
# plot setting
############################################
domain = [-55, 90, 10, 180] 
domain_draw = [-50, 90, 10, 180] 
dlat = 10 
dlon = 30 
#llon, llat = np.meshgrid(lon, lat)
llon=np.copy(lon)
llat=np.copy(lat)
bg_col = '0.6'
cont_col = '1.0'
lev = np.arange(0,1+0.1,0.1)
X,Y,Z=grid(lon,lat,var)

'''
ax=plt.figure(figsize=(10,8)) #12,6)) # (10,8)
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
plt.scatter(llon, llat, var, levels= lev, cmap=plt.cm.afmhot_r) 
cb=plt.colorbar(ticks=lev,shrink=0.9)
#cb.ax1.tick_params(labelsize=14)
plt.title('ppr(ZFL) from cav mtx', fontsize=16, y=1.08)

#plt.savefig(figfile,bbox_inches='tight', format='png', dpi=300)
plt.show(block=False)


plt.figure()
plt.scatter(lon,lat,c=var,vmin=0,vmax=1,cmap=plt.cm.afmhot_r)
plt.colorbar()
plt.show(block=False)
'''


