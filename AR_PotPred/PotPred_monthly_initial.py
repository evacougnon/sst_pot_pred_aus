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

import eac_useful as eac
import eric_oliver as eo

# load useful information
# useful numbers
tau = 6 
# days -- lengths of the chunk
outfile = '/home/ecougnon/ana/PotPred/PP_monthly_SSTa_6m_ZFL'
fname_ = '/home/ecougnon/ana/SSTa_monthly_extremes_Aus'
data_ = np.load(fname_+'.npz')
lon_map = data_['lon_map']
lat_map = data_['lat_map']
SST = data_['SST'].item()
# time period
MinYear = 1982
MaxYear = 2016
NumYears = MaxYear - MinYear+1
NumMonths = NumYears*12
# define the start/end indexes for each chunk
str_id = range(0,NumMonths-1,tau)
end_id = range(tau-1,NumMonths+1,tau)
NumChunk = len(end_id)

# allocate memory
PP_TMM = {}
PP_TMX = {}
PP_TMN = {}
PP_Tp10 = {}
PP_Tp90 = {}
keys = ['Var_interC','Var_noise','Var_slow','p']
'''
describe the keys ..
Var_interC -- total inter chunk variance
Var_noise -- variance of the noise over all chunks
Var_slow -- variance ofthe potentially predictable component
p -- Potential inter-chunk predictability
'''
for key in keys:
    PP_TMM[key] = np.empty((len(lat_map),len(lon_map)))
    PP_TMM[key].fill(np.nan)
    PP_TMX[key] = np.empty((len(lat_map),len(lon_map)))
    PP_TMX[key].fill(np.nan)
    PP_TMN[key] = np.empty((len(lat_map),len(lon_map)))
    PP_TMN[key].fill(np.nan)
    PP_Tp10[key] = np.empty((len(lat_map),len(lon_map)))
    PP_Tp10[key].fill(np.nan)
    PP_Tp90[key] = np.empty((len(lat_map),len(lon_map)))
    PP_Tp90[key].fill(np.nan)

icnt_id = 0
for i in range(0,len(lon_map)):
    t_lon = time.time()
    print(icnt_id+1, 'of', len(lon_map))
    
    jcnt_id = 0
    for j in range(0,len(lat_map)):
# checking that the time series is not empty -- detrend does not
# handle nans... assume that if nans in tmm, same for all the others
        if ((np.isnan(np.min(SST['TMM'][j,i,:])) == False) & \
            (np.isnan(np.max(SST['TMM'][j,i,:])) == False)):
# detrend data -- using the linear least squares fit
            dsst_tmm = signal.detrend(SST['TMM'][j,i,:],type='linear')
            dsst_tmx = signal.detrend(SST['TMX'][j,i,:], type='linear')
            dsst_tmn = signal.detrend(SST['TMN'][j,i,:], type='linear')
            dsst_tp90 = signal.detrend(SST['Tp90'][j,i,:], type='linear')
            dsst_tp10 = signal.detrend(SST['Tp10'][j,i,:], type='linear')
#            '''
## deseasonned
#            dsea_sst_tmm_, season, beta = eo.deseason_harmonic(dsst_tmm,4,12)
#            dsea_sst_tmm = np.squeeze(np.array(dsea_sst_tmm_))
#            dsea_sst_tmx_, season, beta = eo.deseason_harmonic(dsst_tmx,4,12)
#            dsea_sst_tmx = np.squeeze(np.array(dsea_sst_tmx_))
#            dsea_sst_tmn_, season, beta = eo.deseason_harmonic(dsst_tmn,4,12)
#            dsea_sst_tmn = np.squeeze(np.array(dsea_sst_tmn_))
#            dsea_sst_tp90_, season, beta = eo.deseason_harmonic(dsst_tp90,4,12)
#            dsea_sst_tp90 = np.squeeze(np.array(dsea_sst_tp90_))
#            dsea_sst_tp10_, season, beta = eo.deseason_harmonic(dsst_tp10,4,12)
#            dsea_sst_tp10 = np.squeeze(np.array(dsea_sst_tp10_))
#            ''' 

# apply PotPred_ZhengFred function
            [var1, var2, var3, var4] = eac.PotPred_ZhengFred(dsst_tmm, \
                                                             tau,NumChunk, \
                                                             str_id,end_id)

            PP_TMM['Var_interC'][jcnt_id,icnt_id] = var1
            PP_TMM['Var_noise'][jcnt_id,icnt_id] = var2
            PP_TMM['Var_slow'][jcnt_id,icnt_id] = var3
            PP_TMM['p'][jcnt_id,icnt_id] = var4

            [var1, var2, var3, var4] = eac.PotPred_ZhengFred(dsst_tmx, \
                                                             tau,NumChunk, \
                                                             str_id,end_id)
        
            PP_TMX['Var_interC'][jcnt_id,icnt_id] = var1
            PP_TMX['Var_noise'][jcnt_id,icnt_id] = var2
            PP_TMX['Var_slow'][jcnt_id,icnt_id] = var3
            PP_TMX['p'][jcnt_id,icnt_id] = var4

            [var1, var2, var3, var4] = eac.PotPred_ZhengFred(dsst_tmn, \
                                                             tau,NumChunk, \
                                                             str_id,end_id)

            PP_TMN['Var_interC'][jcnt_id,icnt_id] = var1
            PP_TMN['Var_noise'][jcnt_id,icnt_id] = var2
            PP_TMN['Var_slow'][jcnt_id,icnt_id] = var3
            PP_TMN['p'][jcnt_id,icnt_id] = var4

            [var1, var2, var3, var4] = eac.PotPred_ZhengFred(dsst_tp90, \
                                                             tau,NumChunk, \
                                                             str_id,end_id)
    
            PP_Tp90['Var_interC'][jcnt_id,icnt_id] = var1
            PP_Tp90['Var_noise'][jcnt_id,icnt_id] = var2
            PP_Tp90['Var_slow'][jcnt_id,icnt_id] = var3
            PP_Tp90['p'][jcnt_id,icnt_id] = var4

            [var1, var2, var3, var4] = eac.PotPred_ZhengFred(dsst_tp10, \
                                                             tau,NumChunk, \
                                                             str_id,end_id)

            PP_Tp10['Var_interC'][jcnt_id,icnt_id] = var1
            PP_Tp10['Var_noise'][jcnt_id,icnt_id] = var2
            PP_Tp10['Var_slow'][jcnt_id,icnt_id] = var3
            PP_Tp10['p'][jcnt_id,icnt_id] = var4

        jcnt_id = jcnt_id + 1

    elapsed_lon = time.time() - t_lon
    print('elapsed time for each lon:', elapsed_lon)
    icnt_id = icnt_id + 1

# saving data
np.savez(outfile, lat_map=lat_map, lon_map=lon_map, PP_TMM=PP_TMM, \
         PP_TMX=PP_TMX, PP_TMN=PP_TMN, PP_Tp90=PP_Tp90, PP_Tp10=PP_Tp10)
print('data saved')



