'''
	generate the Ningaloo Nino Index from the monthly
	SSTa from the NOAA OISST data set. Anomalie calculated 
	from the 1982-2016 time series

	Author: Eva C
	Created: Nov 2017
'''
# libriries
import numpy as np
from scipy import io # load matlab file
import pandas as pd

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_useful as eac

from matplotlib import pyplot as plt
import matplotlib

#figfile = '/home/ecougnon/ana/EAC_trp.png'
#outfile = '/home/ecougnon/data/EAC_trp'

# load data
matobj = io.loadmat('/home/ecougnon/data/TRANS_EAC.mat')
trp_eac_mat = matobj['transEAC']
date_eac = matobj['dateEAC']

# time vector
MinYear = 1992
MaxYear = 2005
tim_vec = pd.date_range('1992-10-01','2005-12-31',name='time',freq='M')

# allocate memory
trp_eac = np.empty(len(tim_vec))
trp_eac.fill(np.nan)
# check how many sample do we have per month
count_sample = np.empty(len(tim_vec))
count_sample.fill(np.nan)

# make monthly EAC trp mean
k=0
for ii in range(MinYear, MaxYear+1):
    if ii==MinYear:
        for mm in range(10,12+1):
            tmp = np.nanmean(trp_eac_mat[0,np.where((date_eac[:,0] == MinYear) & \
                                                    (date_eac[:,1] == mm))])
            tmp_cnt = np.count_nonzero(trp_eac_mat[0,np.where((date_eac[:,0] == \
                                                               MinYear) & \
                                                              (date_eac[:,1] == mm))])
            trp_eac[k] = tmp
            count_sample[k] = tmp_cnt
            k = k+1
    elif ii>MinYear:
        for mm in range(1, 12+1):
            tmp = np.nanmean(trp_eac_mat[0,np.where((date_eac[:,0] == ii) & \
                                                    (date_eac[:,1] == mm))])
            tmp_cnt = np.count_nonzero(trp_eac_mat[0,np.where((date_eac[:,0] == ii) & \
                                                              (date_eac[:,1] == mm))])
            trp_eac[k] = tmp
            count_sample[k] = tmp_cnt
            k = k+1

#np.savez(outfile, trp_eac=trp_eac, count_sample=count_sample, tim_vec=tim_vec)


plt.figure(figsize=(13,7))
plt.plot(tim_vec, -trp_eac)
plt.plot(tim_vec[(12*2/2):-(12*2/2)+1], eac.moving_average(-trp_eac, 12*2))
plt.grid()
plt.title('Monthly southward EAC transport from Ridgway et al 2008 -- similar to their Fig 13 and Holbrook et al 2011 Fig 2')
plt.ylabel('Transport (Sv)')
plt.ylim([-10, 30])
plt.xlim(['1992-01-01', tim_vec[-1]])
plt.legend(['SSTa', '2 year moving average'])
 
#plt.savefig(figfile, bbox_inches='tight',dpi=300)
plt.show()


'''
# plotting -- maps
domain = [-35, 100, -20, 120] #[-55, 90, 10, 180] 
domain_draw = [-35, 100, -20, 120] #[[-50, 90, 10, 180] 
dlat = 2 #10
dlon = 2 #30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'

ax=plt.figure()
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
plt.contourf(lonproj, latproj, np.nanmean(SST, axis=0), levels=np.arange(-5,5+1,1), \
             cmap=plt.cm.YlOrBr)
'''











