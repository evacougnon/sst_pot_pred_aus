'''

    Compute potential predictability on time series
    of each point following von Storch and Zwiers's book
    
    TESTING PHASE!!!

    run time: about 45 minutes for the region around Oz using
                   35 chunks
              likely to be just under 2 hours for 140 chunks
                   (3 months chunk length)    
              under 1.5 hour for the 6 months chunk length
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
import eric_oliver as eo

# load useful information
pathroot = '/home/ecoliver/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
matobj_tmp = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.0001.mat')
lat_tmp = matobj_tmp['lat']
lon_tmp = matobj_tmp['lon']
# useful numbers
res = 0.25 # resolution of the data -- quarter degree
tau = monthrange(1982,1)[1]+monthrange(1982,2)[1]+monthrange(1982,3)[1] 
# days -- lengths of the chunk
outfile = '/home/ecougnon/ana/AR1/PP_JFM_dsea'
sample = 1 # to use in eac.def_chunk if 0 chunk follow
           # each other if 1 every year sampling (season -- default winter)
# time period
MinYear = 1982
MaxYear = 2016
NumYears = MaxYear - MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
dtime = [datetime.date(MinYear,1,1) + datetime.timedelta(days=i) \
         for i in range(NumDays)]
# remove the 29 of Feb
idx_leap=np.zeros(MaxNumLeapYear)
k=0
for t in range(0,len(dtime)):
    if ((dtime[t].month==2) &  (dtime[t].day==29)):
        idx_leap[k]=t
        k = k + 1
dtime=np.delete(dtime,idx_leap)    
# define the start/end indexes for each chunk
NumChunk = NumYears #int(len(dtime)/tau) 
# number of 91 day-chunk within the time series 
# NumCunck = NumYears #if looking at a season every year
[str_id, end_id]=eac.def_chunk(dtime, tau, NumChunk, sample, \
                               str_month=1, str_day=1, \
                               end_month=3, end_day=31)

# define region of study 
lat_min = -55 #-17 #-55
lat_max = 10.25 #0 #10.25
lon_min = 90 #155 #90
lon_max = 180.25 #165 #180.25
# find the index of the lat/lon min/max
lat_min_id=np.nonzero((lat_tmp>(lat_min-res)) & (lat_tmp<(lat_min)))
lat_max_id=np.nonzero((lat_tmp>(lat_max-res)) & (lat_tmp<(lat_max)))
lon_min_id=np.nonzero((lon_tmp>(lon_min-res)) & (lon_tmp<(lon_min)))
lon_max_id=np.nonzero((lon_tmp>(lon_max-res)) & (lon_tmp<(lon_max)))
# get the matrix of indexes for the lat/lon
lat_map = lat_tmp[lat_min_id[0][0]:lat_max_id[0][0]+1]
lat_id = range(lat_min_id[0][0],lat_max_id[0][0]+1,1)
lon_map = lon_tmp[lon_min_id[0][0]:lon_max_id[0][0]+1]
lon_id = range(lon_min_id[0][0],lon_max_id[0][0]+1,1)
# allocate memory
PP = {}
keys = ['Var_interC','Var_F','S_tau','F90','F95']
'''
Var_interC -- inter-chunk variability
Var_F -- inter chunk variance due to the fast variaying component
S_tau -- relative importance of the potentially predictable
         component
F90
F95
'''
for key in keys:
    PP[key] = np.empty((len(lat_map),len(lon_map)))
    PP[key].fill(np.nan)

icnt_id = 0
for i in lon_id:
    t_lon = time.time()
    print(icnt_id+1, 'of', len(lon_id))
# load time series from the given point
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' \
                        + str(i+1).zfill(4) + '.mat')
                        # (i+1) due to difference between matlab
                        # and python indexing
    jcnt_id = 0
    for j in lat_id:
#        t_lat = time.time()
        sst_ts = matobj['sst_ts'][j][:]
# remove the 29 of Feb from the time series
        sst_ts_ = np.delete(sst_ts, idx_leap)
# checking if continent or not -- does not handle empty matrix
        if ((np.isnan(np.min(sst_ts_)) == False) & \
            (np.isnan(np.max(sst_ts_)) == False)):

# detrend the data -- using the linear least squares fit
            dsst = signal.detrend(sst_ts_, axis=0, type='linear')
# deseasonned
            dsea_sst, season, beta = eo.deseason_harmonic(dsst,6,365)
# allocate memory
## time series including only the chunks
#            dsst_chunk=np.empty(tau*NumChunk)
#            dsst_chunk.fill(np.nan)
## filtered chunk time series
#            dsst_box=np.empty(len(dsst_chunk))
#            dsst_box.fill(np.nan)
# chunk mean allocation
            t_j=np.empty(NumChunk) 
            t_j.fill(np.nan)
# periodogram of each chunk
            pxx_j = np.empty([NumChunk,int(tau/2)+1]) 
            pxx_j.fill(np.nan)# power spectrum of each chunk
            f_j =np.empty([NumChunk,int(tau/2)+1]) 
            f_j.fill(np.nan)# corresponding frequencies
#            pxx_box = np.empty([NumChunk,int(tau/2)+1])
#            pxx_box.fill(np.nan)# power spectrum of each chunk
#            f_box =np.empty([NumChunk,int(tau/2)+1])
#            f_box.fill(np.nan)# corresponding frequencies
#            t_spec = time.time()
            kk = 0
            for c in range(0,NumChunk):
# chunk mean
                t_j[c] = np.nanmean(dsea_sst[int(str_id[c]):int(end_id[c])+1])
# spectrum of each chunck for the fast varying component
                f_j[c,:],pxx_j[c,:]=signal.periodogram(dsea_sst[int(str_id[c]) \
                                                                :int(end_id[c])+1]) 
#            elapsed_spec = time.time() - t_spec
#            print('elapsed time for each spectrum calc:', elapsed_spec)
## time series of the chunks                
#                dsst_chunk[kk:kk+tau] = dsst[int(str_id[c]) \
#                                             :int(end_id[c])+1]
                kk = kk + tau
# mean of all the chunks
            t_jall = np.nanmean(t_j)

# np.var(t_j) # inter chunck varinance is different as it uses
# 1/(n-1) instead of 1/n in the variance calculation
# commun practice -- sample variance
# inter chunk variance (sample variance)
            tmp=np.zeros(NumChunk)
            for c in range(0,NumChunk):
                tmp[c] = (t_j[c]-t_jall)**2
            PP['Var_interC'][jcnt_id,icnt_id] = 1/(NumChunk-1) * np.nansum(tmp)
### running mean filter of the combined chunk time series -- 
### with a window of the length of the chunk -- keep the high frequencies
#            window = tau
#            weights = np.repeat(1.0,window)/window
#            dsst_box = np.convolve(dsst_chunk, weights,'valid')
### power spectral density of the chunck time series low pass filtered
#            nn = 0
#            for c in range(0,NumChunk-1):
#                f_box[c,:], pxx_box[c,:] = signal.periodogram(dsst_box[nn \
#                                                                       :nn+tau])
#                nn = nn + tau

# inter chunk variance of the fast varying process
# average periodograms in time -- estimate the spectrum for all the chunks
            pxx_jall=np.mean(pxx_j,axis=0)

            idx = np.argmin(np.abs(f_j[0][:] - (1/tau)))
# variance of the fast varying component
            PP['Var_F'][jcnt_id,icnt_id] = (1/tau) * pxx_jall[idx]
#            PP['Var_F2'][jcnt_id,icnt_id] = np.nanmean(pxx_box[:][idx]/tau)

        jcnt_id = jcnt_id + 1
#        elapsed_lat = time.time() - t_lat
#        print('elapsed time for each lat:', elapsed_lat)

    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt_id = icnt_id + 1

# variance ratio
PP['S_tau'][:,:] = PP['Var_interC'][:,:] / PP['Var_F'][:,:]
#PP['S_tau2'][:,:] = PP['Var_interC'][:,:] / PP['Var_F2'][:,:]
##
# testing the null hypothesis to find the significance
PP['F90'] = stats.f.ppf(0.90,tau-1,2*NumChunk)
PP['F95'] = stats.f.ppf(0.95,tau-1,2*NumChunk)
PP['F10'] = stats.f.ppf(0.10,tau-1,2*NumChunk)
PP['F05'] = stats.f.ppf(0.05,tau-1,2*NumChunk)


# saving data
np.savez(outfile, lat_map=lat_map, lon_map=lon_map, PP=PP)
print('data saved')

'''
# plotting test
# plot setting
domain = [-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180]
dlat = 10
dlon = 30
llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'

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
plt.contourf(lonproj, latproj, PP['S_tau'], levels=np.arange(0,6+0.1,0.1), \
             cmap=plt.cm.viridis)
plt.colorbar()
plt.contour(lonproj, latproj, PP['S_tau'], levels=[PP['F95']], colors='1')
plt.title('Potential Pred. -- JFM -- Zwiers, von S. method')
# -- with F95= %d and F90 = %d' %(PP['F95'],PP['F90']))
plt.show()
'''


