# function to calc the time series of a climate index

import numpy as np
import scipy as sp
from scipy import stats
from scipy import signal
import pandas as pd
import time as time
import datetime
import xarray as xr
from glob import glob

import sys
sys.path.insert(0,'/home/ecougnon/scripts/libraries/')
import eac_useful as eac

def calc_month_index(index, MinYear, MaxYear) :
    '''
    return the monthly index of the specified climte index 
    with the trend removed and standardise index 
    (monthly - mean)/std
    uses 'commun' index

    INPUT:
    index -- name of the index wanted, can be:
        'MEI'
        'NINO34'
        'DMI' for the IOD
        'SAM'
        'NNI' for Ningaloo Nino index from Ming Feng paper -- check date!
    MinYear
    MaxYear

    OUTPUT:
    ind_monthly -- monthly time series of the index detrended
    ind_monthly_std -- same than above but standardised   

    '''

    if ((index=='MEI') | (index=='NINO34') | (index=='DMI') | (index=='SAM')):

        header = '/home/data/index/'
        NumYears = MaxYear-MinYear+1
    
        ind_monthly = np.empty(NumYears*12)
        ind_monthly_std = np.empty(NumYears*12)
        k=0

        if index=='MEI':
          file_ind = np.genfromtxt(header + 'enso/MEI_index.txt', \
                                   skip_header=10, skip_footer = 30, \
                                   delimiter='\t')
        elif index=='NINO34':
            file_ind = np.genfromtxt(header + 'enso/nino34.txt', \
                                     skip_header=1, skip_footer = 6)
        elif index=='DMI':
            file_ind = np.genfromtxt(header + 'DMI_IOD.txt', \
                                     skip_header=1, skip_footer = 4)
        elif index=='SAM':
            file_ind = np.genfromtxt(header + 'sam/marshall2003.txt', \
                                     skip_header=2, skip_footer = 3)

        str_ind = np.nonzero((file_ind[:,0]>(MinYear-1)) \
                             & (file_ind[:,0]<(MinYear+1)))
        for yy in np.arange(str_ind[0][0],len(file_ind[:,0])):
            for mm in np.arange(1,12+1):
                ind_monthly[k] = file_ind[yy,mm]
                k = k + 1
        ind_monthly = signal.detrend(ind_monthly)
        ind_monthly_std = (ind_monthly-np.nanmean(ind_monthly)) \
                          /np.nanstd(ind_monthly)

        return ind_monthly, ind_monthly_std

    elif index=='NNI':
        data_NNI = np.load('/home/ecougnon/data/NingalooNino_index.npz')
        ind_monthly = data_NNI['NingN_Ind']
        ind_monthly = signal.detrend(ind_monthly)
        ind_monthly_std = (ind_monthly-np.nanmean(ind_monthly)) \
                          /np.nanstd(ind_monthly)

        return ind_monthly, ind_monthly_std

    else:
        return print('WARNING: the index is not defined in the function, update the function or change the index. Also note that EAC transport index is in another function: calc_month_eac')



def calc_month_eac(tim_vec):
    '''
    calc the monthly EAC transport from ??? check which paper 

    INPUT:
    time_vec -- timeseries of the variable to compare with EAC trp index
    
    OUTPUT:
    eac_str, eac_end -- the corresponding stat and end index of time_vec
                        that matches with the EAC time series
    ind_monthly -- monthly time series
    '''


    data_eac= np.load('/home/ecougnon/data/EAC_trp.npz')
    eac_time = data_eac['tim_vec']
    eac_monthly = data_eac['trp_eac']
    eac_str = int(np.array(np.where(tim_vec == eac_time[0])))
    eac_end = int(np.array(np.where(tim_vec == eac_time[-1])))
    eac_monthly = signal.detrend(eac_monthly)

    return eac_monthly, eac_str, eac_end

def nino34_index_dps():
    '''
    calc NINO3.4 index with the DPS
    using enkf-9 (reanalysis)

    output the SSTa averaged for the NINO 3.4 region
    and the phase: +ve(+1), neutral 90) and -ve (-1) following:
    
    An El Nino or La Nina event is identified if the 5-month 
    running-average of the NINO3.4 index exceeds +0.4$^\circ$C for 
    El Nino or -0.4$^\circ$C for La Nina for at least 6 consecutive months.
    '''

    header_out = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
    header = '/home/ecougnon/data/DPS/reanalysis/ETKF/'
    gname = header + 'ocean_daily_SST_enkf9_mem001_2002.nc'
    outfile = header_out + 'nino34_enkf9_mem001_20032017.nc'
# define the region 
    lat_min = -5
    lat_max = 5
# !!! LONGITUDE in the model!!! from -280 to 80 deg E
    lon_min = 190 - 360
    lon_max = 240 - 360
    yt_ocean = xr.open_dataset(gname)['yt_ocean']
    yt_ocean = yt_ocean.sel(yt_ocean=slice(lat_min,lat_max))
    xt_ocean = xr.open_dataset(gname)['xt_ocean'] #!!! from -280 to 80 deg E
    xt_ocean = xt_ocean.sel(xt_ocean=slice(lon_min,lon_max))
# monthly mean for the Australian region
    SST_m = eac.read_netcdfs('/home/ecougnon/data/DPS/reanalysis/ETKF/ocean_daily_SST_enkf9_mem001_20??-20??.nc', \
                             dim='time', lat=yt_ocean, lon=xt_ocean, \
                             transform_func=lambda ds: ds.resample('1MS', \
                                                                   dim='time', \
                                                                   how='mean'))
    SST_m = SST_m.squeeze('st_ocean')
    SST_m = SST_m.sel(time=slice('2003-01-01','2017-11-30'))
    SST_m = SST_m.mean(dim=('yt_ocean','xt_ocean'))
# calc monthly climatology
    clim_m = SST_m.groupby('time.month').mean('time')
    ssta_m = SST_m.groupby('time.month') - clim_m
#    dssta_m = np.apply_along_axis(signal.detrend,0,ssta_m['temp'], type='linear')
#    nino34_id = np.mean(dssta_m,axis=(1,2))
    nino34_id = np.array(ssta_m['temp'])
# check +ve and -ve phases (-1 = -ve; 0 = neutral; +1 = +ve)
    nino34_id_smoothed = pd.DataFrame(nino34_id). \
                            rolling(5,center=True,win_type='boxcar',min_periods=5). \
                            sum()/5
    nino34_id_smoothed = np.squeeze(np.array(nino34_id_smoothed))
    nino34_ph = np.zeros(nino34_id_smoothed.shape)
    nino34_ph[nino34_id_smoothed>=0.4] = 1
    nino34_ph[nino34_id_smoothed<=-0.4] = -1
    
    test = nino34_ph.copy()

    tt=0
    while tt < len(nino34_ph):
        if nino34_ph[tt] == 0:
            tt = tt+1
        elif ((nino34_ph[tt] == 1) & (np.sum(nino34_ph[tt:tt+6]) < 6)):
            nino34_ph[tt:tt+6] = 0
            tt = tt+6
        elif ((nino34_ph[tt] == 1) & (np.sum(nino34_ph[tt:tt+6]) == 6)):
            tt = tt+6
            while nino34_ph[tt] == 1:
                tt = tt+1
        elif ((nino34_ph[tt] == -1) & (np.sum(nino34_ph[tt:tt+6]) > -6)):
            nino34_ph[tt:tt+6] = 0
            tt = tt+6
        elif ((nino34_ph[tt] == -1) & (np.sum(nino34_ph[tt:tt+6]) == -6)):
            tt = tt+6
            while nino34_ph[tt] == -1:
                tt = tt+1


    return nino34_id,nino34_ph


