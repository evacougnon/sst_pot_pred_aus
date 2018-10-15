'''
    plot the power spectrum for EOF analysis
    giving the time band variance

    Created: Oct 2018
    Author: Eva C.
'''
# import libraries
import numpy as np
import xarray as xr
import pandas as pd
from scipy import signal
import scipy.stats.mstats as st # highly similar to scipy.stats 
                                # but adapted for masked arrays

import sys
sys.path.insert(0,'../libraries/')
import eac_useful as eac

import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean

# Load data
fname = '../../ana/PotPred/EOF/eof_Aus_daily_trend_1deg_12modes_monthly.nc'
figfile_psd ='/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/PCsm_psd_5-8modes_SSTatrend.png'
figfile_psdall ='/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/PCsm_psdall_5-8modes_SSTatrend_zoom.png'

lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
time_ts = xr.open_dataset(fname)['time']
PCs = xr.open_dataset(fname)['PCs']


###########################################
## spectral analysis on the PCs
###########################################
fs=1
f_dblside, pxx_dblside = np.apply_along_axis(signal.periodogram,0, \
                                             PCs.transpose('time','modes'),fs, \
                                             window ='hanning',detrend=False, \
                                             return_onesided=False, \
                                             scaling='density')
window = 3
M = PCs.shape[0]
f_dblside_ = np.nan*np.zeros((len(f_dblside[:,0]),M))
pxx_dblside_ = np.nan*np.zeros((len(f_dblside[:,0]),M))
pxx_dblside_smooth = np.nan*np.zeros((len(f_dblside[:,0]),M))
for ii in range(0,M):
    id_sort = np.argsort(f_dblside[:,ii],axis=0)
    f_dblside_[:,ii] = f_dblside[id_sort,ii]
    pxx_dblside_[:,ii] = pxx_dblside[id_sort,ii]*2
    pxx_dblside_[int(len(pxx_dblside_)/2),ii] = \
                pxx_dblside_[int(len(pxx_dblside_)/2),ii]/2
    pxx_dblside_smooth[int(window/2):-int(window/2),ii]=  \
                       eac.movingaverage(pxx_dblside_[:,ii],window)

InBand0_p = np.nan*np.zeros(M) # band between <1 years
InBand1_p = np.nan*np.zeros(M) # band between 1-2 years or 1-3
InBand2_p = np.nan*np.zeros(M) # band between 2-4 years
InBand4_p = np.nan*np.zeros(M) # band between 4-7 years or 3-7
InBand7_p = np.nan*np.zeros(M) # band between 7-10 years
InBand10_p = np.nan*np.zeros(M) # band between >10 years
InBand_all_p = np.nan*np.zeros(M) # band between all
pxx = pxx_dblside_smooth[int(window/2):-int(window/2),:] #pxx_pwe #pxx_per
f = f_dblside_[int(window/2):-int(window/2),:] #f_pwe #f_per
max_f = np.max(np.max(f,axis=0))
for m in range(0,M):
    InBand0_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0.08,max_f)
    InBand1_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0.0278,0.08)
#0.042,0.08)
    InBand2_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0.021,0.042)
    InBand4_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0.012,0.0278)
#0.012,0.021)
    InBand7_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0.0083,0.012)
    InBand10_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0,0.0083)
    InBand_all_p[m] = eac.bandpower_Pxx(pxx[:,m],f[:,m],0,max_f)

################################################
# normalised by the total variance of each mode
################################################
InBand0_p_norm = 100 * (InBand0_p / InBand_all_p)
InBand1_p_norm = 100 * (InBand1_p / InBand_all_p)
InBand2_p_norm = 100 * (InBand2_p / InBand_all_p)
InBand4_p_norm = 100 * (InBand4_p / InBand_all_p)
InBand7_p_norm = 100 * (InBand7_p / InBand_all_p)
InBand10_p_norm = 100 * (InBand10_p / InBand_all_p)


'''
###############
# plot all psd
################
hmode = 1
plt.figure(figsize=(17,4))
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-hmode], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-hmode],'*-')
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+1)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+1)],'*-')
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+2)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+2)],'*-')
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+3)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+3)],'*-')
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+4)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+4)],'*-')
plt.semilogy(f_dblside_[int(len(pxx_dblside_)/2):,-(hmode+5)], \
             pxx_dblside_smooth[int(len(pxx_dblside_)/2):,-(hmode+5)],'*-')
plt.plot([0.0083, 0.0083],[0, 10**7],'k--') #10 yr
plt.plot([0.012, 0.012],[0, 10**7],'k--') #7 yr
#plt.plot([0.0167, 0.0167],[0, 10**7],'k--') #5yr
plt.plot([0.021, 0.021],[0, 10**7],'k--') # 4yr
#plt.plot([0.0278, 0.0278],[0, 10**7],'k--') # 3yr
plt.plot([0.042, 0.042],[0, 10**7],'k--') # 2yr
plt.plot([0.08, 0.08],[0, 10**7],'k--') # 1yr
plt.plot([0.167, 0.167],[0, 10**7],'k--',alpha=0.3) # 6 months
plt.xlim([0,0.2])
plt.ylim([0,10**5])
plt.grid()
plt.legend(['PC' + str(hmode) + '','PC' + str(hmode+1) + '', \
            'PC' + str(hmode+2) + '','PC' + str(hmode+3) + '', \
            'PC' + str(hmode+4) + '','PC' + str(hmode+5) + ''])
plt.title('PSD for each PCs with 10, 7, 4, 2, 1, 0.5 year marks')
plt.savefig(figfile_psdall, bbox_inches='tight', format='png', dpi=300)
plt.show(block=False)


'''


