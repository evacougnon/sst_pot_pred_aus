'''
    testing PotPred_ts_AR1 with Jiale Lou data
    and thinking a way to write the main code of 
    PotPred_ts_AR1 as a function
'''

# import libraries
import numpy as np
from scipy import signal
from scipy import io
from scipy import stats
import time as time
import datetime
from dateutil.relativedelta import relativedelta
from calendar import monthrange

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import eac_useful as eac


# load data
fname = '/home/ecougnon/data/CCSM4LM_LOU/zg500_mat.mat'
# issues downloading the raw file as some lines are missing value
# opened and saved with matlab 
#         5x5.zg_500_Amon_CCSM4_past1000_r1i1p1_085001-185012.dat'
matobj = io.loadmat(fname)
tmp=matobj['zg500']
# reshape and remove the NaNs
tmp2 = np.reshape(tmp,-1)
zg500_ = np.delete(tmp2,np.where(np.isnan(tmp2)==True))
zg500 = np.reshape(zg500_,(1001,12,37,72))
# use north hemisphere only
zg500_north = zg500[:,:,18:,:]
# interpolation in longitude
# create a masking matrix with ones and zeros
'''
# do  not use the mask for now
mask_mtx=np.ones((72,19)) #,dtype=bool)
mask_mtx[0,-1]=1 # at 90N take a single value
mask_mtx[range(0,72,12),-2]=1 # at 85N take 1 value every 12 lon (6 lon in tot)
mask_mtx[range(0,72,6),-3]=1 # at 80N sample every 30deg of lon
mask_mtx[range(0,72,4),-4]=1 # at 75N sample every 20deg of lon
mask_mtx[range(0,72,3),-5]=1 # at 70N sample every 15deg of lon
mask_mtx[range(0,72,2),-9:-6+1]=1 # 50-65N sample every 10deg of lon
mask_mtx[:,0:-10+1]=1 # 50-65N sample every 10deg of lon
mask_swap=np.swapaxes(mask_mtx,0,1)
mask_vec=np.reshape(mask_swap,19*72)
zg_tmp0=np.swapaxes(zg500_north,0,1)
zg_tmp = np.reshape(zg_tmp0,(19*72,12,1001))
zg500_N_sample = np.empty((19*72,12,1001))
zg500_N_sample.fill(np.nan)
for i in range(0,72*19):
    if (mask_vec[i]==1):
        zg500_N_sample[i,:,:] = zg_tmp[i,:,:]
#zg500_N_sample = np.compress(mask_vec,zg_tmp,axis=0)
# annual mean
zg500_N_mean = np.nanmean(zg500_N_sample[:,:,0:-1],axis=1)
# anomaly
zg_mean=np.nanmean(zg500_N_mean,axis=1)
zg500_N_a=np.zeros((len(zg500_N_mean[:,1]),1000))
for i in range(0,1000):
    zg500_N_a[:,i] = zg500_N_mean[:,i] - zg_mean[:]
'''
# instead of using the mask
# annual mean
zg500_N_mean = np.nanmean(zg500_north[0:-1,:,:,:],axis=1)
# anomaly
zg_mean=np.nanmean(zg500_N_mean,axis=0)
zg500_N_a_=np.zeros((1000,19,72))
for i in range(0,1000):
    zg500_N_a_[i,:,:] = zg500_N_mean[i,:,:] - zg_mean[:,:]
# easier to work on vector
zg500_N_a=np.reshape(zg500_N_a_,(1000,19*72))

# useful value for the potential predictability calc
tau=10 # length of the chunk
NumChunk=100 # total number of chunk
str_id=np.arange(0,1000,10)
end_id=np.arange(9,1000+1,10)
 
# allocate memory for the outputs of the function PotPred
PP = {}
keys = ['Var_interC','Var_noise','Var_slow','p']
'''
describe the keys ..
Var_interC -- total inter chunk variance
Var_noise -- variance of the noise over all chunks
Var_slow -- variance ofthe potentially predictable component
p -- Potential inter-chunk predictability
'''
for key in keys:
    PP[key] = np.empty(19*72)
    PP[key].fill(np.nan)

# for each location
icnt_id = 0
for i in range(0,19*72):
    t_loc = time.time()

# calling the potential pred function
#    [var1, var2, var3, var4] = PotPred_ZhengFred(zg500_N_a[:,i],tau,NumChunk \
#                                                 , str_id,end_id)
    [var1, var2, var3, var4, var5] = eac.PotPred_vStorchZwiers(zg500_N_a[:,i], \
                                                               tau,NumChunk, \
                                                               str_id,end_id)
    PP['Var_interC'][icnt_id] = var1
    PP['Var_noise'][icnt_id] = var2
    PP['p'][icnt_id] = var3
    PP['F90'] = var4
    PP['F95'] = var5

    elapsed_loc = time.time() - t_loc
#    print('elapsed time for each location:', elapsed_loc)
    icnt_id = icnt_id + 1

test=np.reshape(PP['p'],(19,72))
test2=np.reshape(PP['Var_interC'],(19,72))
test3=np.reshape(PP['Var_noise'],(19,72))
#test4=np.reshape(PP['Var_slow'],(19,72))

lon=np.arange(270,-90,-5)
lat=np.arange(0,90+1,5)

# plotting
domain = [0, 0, 90, 360]
domain_draw = [0, 0, 90, 360]
dlat = 20
dlon = 30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'

plt.figure(figsize=(10,8))
plt.clf()
proj = bm.Basemap(projection='npstere',boundinglat=0,lon_0=180,resolution='l')
#proj = bm.Basemap(resolution='c',projection='ortho',lat_0=90.,lon_0=180.)
#, llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#proj.fillcontinents(color=(1,1,1), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawcoastlines()
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[False,False,False,False])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[True,True,True,True])
lonproj, latproj = proj(llon, llat)
plt.contourf(latproj, lonproj,test, levels=np.arange(0.1,0.8+0.1,0.1), \
#plt.contourf(latproj, lonproj,test2, levels=[0,20,50,80,110,140,180], \
             cmap=plt.cm.Oranges)
plt.colorbar()
#plt.title('Potential predictability')
# -- with F95= %d and F90 = %d' %(PP['F95'],PP['F90']))
plt.show()

