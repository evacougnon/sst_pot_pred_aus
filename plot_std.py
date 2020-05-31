# plt std or variance

import numpy as np
import xarray as xr

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm


# Load data

fname = '/v_Munk_Drive/ecougnon/ana/SSTa_daily_Aus.nc'
figfile = '/v_Munk_Drive/ecougnon/figures/SSTa_daily_var.eps'

lat_map = xr.open_dataset(fname)['lat']
lon_map = xr.open_dataset(fname)['lon']
SSTa  = xr.open_dataset(fname)['SSTa']

TMM_std = np.nanvar(SSTa,axis=0)

# plot setting
domain = [-55, 90, 10, 180] #[-80, -180, 85, 180] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon_map, lat_map)
#llat, llon = np.meshgrid(lat_map, lon_map)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure(figsize=(19,7)) #15,5)) #(19,7))
plt.clf()

plt.subplot(1,1,1, facecolor=bg_col)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon_, llat_)
plt.contourf(lonproj, latproj, TMM_std, levels=np.arange(0,1.2+0.1,0.1), \
             cmap=plt.cm.afmhot_r)
cb=plt.colorbar(ticks=np.arange(0,1.2+0.1,0.2),shrink=0.8)
cb.ax.tick_params(labelsize=14)
plt.title('variance daily SSTa', \
          fontsize=14, y=1.02)

plt.savefig(figfile,bbox_inches='tight', format='eps', dpi=300)
plt.show()

