
# Load required modules

import numpy as np
from scipy import io
from scipy import linalg
from scipy import stats
from datetime import date
from netCDF4 import Dataset

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm



fname = '/home/ecougnon/ana/AR1/SST_annual_extremes_Aus'

data = np.load(fname + '.npz')
lon_map = data['lon_map']
lat_map = data['lat_map']
years = data['years']
SST = data['SST'].item()

# Re-map to run 20E to 380E
i_20E = np.where(lon_map>20)[0][0]
lon_map = np.append(lon_map[i_20E:], lon_map[:i_20E]+360)
key = 'TMN' 
SST[key] = np.append(SST[key][:,i_20E:,:], SST[key][:,:i_20E,:], axis=1)

SST_mean = np.nanmean(SST[key],axis=2)

# plot setting
domain = [-18, 150, -15, 158] #[-50, 90, 0, 180]
domain_draw = [-18, 150, -15, 158] # [-50, 90, 0, 180]
dlat = 1 #20
dlon = 2 #30
llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'

# plotting
plt.figure(figsize=(10,8))
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
		  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
		   labels=[True,False,False,False])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
		   labels=[False,False,False,True])
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, SST_mean, levels=np.arange(-2,15+0.5,0.5), \
	     cmap=plt.cm.viridis)
H = plt.colorbar()
plt.title('Mean annual min SST')
plt.scatter(2,2,s=10,c='k', alpha=1)



