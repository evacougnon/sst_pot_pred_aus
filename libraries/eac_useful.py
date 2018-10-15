# function to find the start and end index of a chunck
# for the potential predictability script

import numpy as np
import scipy as sp
from scipy import stats
from scipy import signal
import time as time
import datetime
import xarray as xr
from glob import glob
from geopy.distance import vincenty # distance lat/lon

def find_nearest(array,value):
    '''
    find the nearest value from the array and return the index within the 
    array and the value
    '''
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def dxdy(lat,lon) : 
    '''
    Uses geopy.distance.vincenty to calculate the distance between each
    grid cell -- the default ellipsoid is WGS84
    Made for oservational data where the data point is the centre of the cell

    INPUT:
    lat, lon -- vectors

    OUTPUT:
    dx, dy -- distance for each cell edge in m
    '''

    dx_tmp = np.NaN*np.zeros((len(lat)-1, \
                              len(lon)-1))
    dy_tmp = np.NaN*np.zeros((len(lat)-1, \
                              len(lon)-1))
    for lt in np.arange(0,len(lat)-1):
        for ln in np.arange(0,len(lon)-1):
            lat_tmp = lat[lt]
            lon_tmp = lon[ln]
            coord1 = (lat_tmp,lon[ln])
            coord2 = (lat_tmp,lon[ln+1])
            coord1_ = (lat[lt],lon_tmp)
            coord2_ = (lat[lt+1],lon_tmp)

            dx_tmp[lt,ln] = (vincenty(coord1,coord2).km)/2
            dy_tmp[lt,ln] = (vincenty(coord1_,coord2_).km)/2

    dx = np.NaN*np.zeros((len(lat)-2, \
                          len(lon)-2))
    dy = np.NaN*np.zeros((len(lat)-2, \
                          len(lon)-2))
    for lt in np.arange(0,len(lat)-2):
        for ln in np.arange(0,len(lon)-2):
            dx[lt,ln] = (dx_tmp[lt,ln] + dx_tmp[lt,ln+1]) * 1000
            dy[lt,ln] = (dy_tmp[lt,ln] + dy_tmp[lt+1,ln]) * 1000

    return dx, dy
'''
def moving_average(ts, w) :
    ''
    Low pass filter type: box-car or moving averaged
    INPUT:
    ts -- time series
    w -- window size of the box for averaged
    OUTPUT:
    return the low pass filtered time series with a shorter size:
    len(ts_output) = ts[(w/2):-(w/2)+1]

    NOTES: gives the same results than using 
    np.convolve(ts, weights, 'valid') with weights = np.repeat(1.0,w)/w
    '' 
    ret = np.cumsum(ts, dtype=float)
    ret[w:] = ret[w:] - ret[:-w]

    return ret[w - 1:] / w
'''
def movingaverage(values, window):
    '''
    Input:
        values: time series
        window: number of points to include in the window

    Output:
        MA: the time series resulting of the moving average with 
            a shorter length
    '''
    weights = np.repeat(1.0, window)/window
    MA = np.convolve(values, weights, 'valid')
    return MA

def get_autocorr_coef(ts1,ts2):
    '''
    When using the np.corrcoef, returns the correlation matrix of the 
    variables. In the case of the auto correlation, only one value is 
    wanted (alpha in this function)
    INPUT:
    ts1 and ts2 -- same as np.corrcoef -- the 2 time series to correlate
                for auto-corr of lag one ts[1:] and ts[0:-1]
    OUTPUT:
    alpha -- the auto correlation lag 1 coefficient
    '''
# check for nans
    valid = ~np.isnan(ts1)
# checking that the matrix is not entirely nans
    if (valid.any()==True):
# conpute the auto-correlation coef
        corr_coef_mtx = np.corrcoef(ts1,ts2)
        alpha = corr_coef_mtx[0][1]
# if the time series is empty (only nans), return np.nan
    elif (valid.all()==False):
        alpha = np.nan

    return alpha

##################################################
# in band varance analysis
##################################################
def bandpower(x, fs, fmin, fmax):
    '''
    Taken (and adapted -- to not miss indexes) from the internet: 
    https://stackoverflow.com/questions/44547669/python-numpy-equivalent-of-bandpower-from-matlab
    Based on the matlab function (bandpower)

    INPUT:
    x - time series
    fs - sample rate to return the power in a specified frequency band
    fmin - lower band of the frequency range
    fmax - upper band of the frequency range

    Return the average power in the frequency range
    '''
    f, Pxx = signal.periodogram(x, fs=fs)
    ind_min = sp.argmax(f > fmin) - 1
    ind_max = sp.argmax(f > fmax) - 1
    if ind_max != 0:
        var = sp.trapz(Pxx[ind_min: ind_max+1], f[ind_min: ind_max+1])
    elif ind_max == 0:
        var = sp.trapz(Pxx[ind_min:], f[ind_min:])
    return var

def bandpower_Pxx(Pxx, f, fmin, fmax):
    ''' 
    Adapted from (see below) to read directly the power spectrum

    Taken (and adapted -- to not miss indexes) from the internet: 
    https://stackoverflow.com/questions/44547669/python-numpy-equivalent-of-bandpower-from-matlab
    Based on the matlab function (bandpower)

    INPUT:
    Pxx - power spectrum
    f - fequencies from the power spectrum
    fmin - lower band of the frequency range
    fmax - upper band of the frequency range

    Return the average power in the frequency range by integrating 
        the power spectral density (PSD) estimate, pxx (usingthe  trapezoidal rule)
    '''
    ind_min = sp.argmax(f > fmin) - 1
    ind_max = sp.argmax(f > fmax) - 1
    if ind_max != -1: 
        var = sp.trapz(Pxx[ind_min: ind_max+1], f[ind_min: ind_max+1])
    elif ind_max == -1:
        var = sp.trapz(Pxx[ind_min:], f[ind_min:])
    return var

####################################################
##### reading and combining netcdf files with the 
##### xarray library and eventually apply a simple
##### transformation along the files (mean for instance)
##### Added in Nov 2017
#####################################################
def read_netcdfs(files, dim, lat, lon, transform_func=None):
    '''
    reading, slicing in lat/lon and combining netcdf files with the 
    xarray library and eventually apply a simple
    transformation along the files (mean for instance)
    Added in Nov 2017

    !!!!!!!!!!!!!!!!!!! 
    for now (Nov 2017) work only for the a latitude variable 
    named yt_ocean and a longitude variable named xt_ocean as 
    in the simulations ran by Terry O'Kane group
    !!!!!!!!!!!!!!!!!!!!!!!    

    INPUT:
    files    -- path to files
    dim      -- name the dimention that the data will be combined
		along
    lat      -- latitude slice for subsetting the data (yt_ocean)
    lon      -- longitude slice for subsetting the data (xt_ocean)

    OUTPUT
    combined -- the combined dataset

    Here we suppose we only care about the combined mean of each file;
    you might also use indexing operations like .sel to subset datasets
    example below:
    combined = read_netcdfs('/all/my/files/*.nc', dim='time',
                            transform_func=lambda ds: ds.mean())

    copied from the xarray website:
       http://xarray.pydata.org/en/stable/io.html

    '''
    def process_one_path(path):
        # use a context manager, to ensure the file gets closed after use
        with xr.open_dataset(path) as ds:
            # slicing/subsampling the data
            ds = ds.sel(xt_ocean=lon, yt_ocean=lat)
            # transform_func should do some sort of selection or
            # aggregation
            if transform_func is not None:
                ds = transform_func(ds)
            # load all data from the transformed dataset, to ensure we can
            # use it after closing each original file
            ds.load()
            return ds

    paths = sorted(glob(files))
    datasets = [process_one_path(p) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined


####################################################
##### potential predictability decomposition related
##### ~~ 2017 (Apr-Jul)
#####################################################
def def_chunk(dtime, tau, NumChunk, sample, MinYear=1982, \
              MaxYear=2016, str_month=6, str_day=1, \
              end_month=8, end_day=31):

    '''
    Input:
    dtime -- time vector of the time series in datetime.date format
    tau -- length of a chunk (number of points within a chunck)
    NumChunk -- number of chunk for the time series
    sample -- 0 each chunk follow each other
              1 changes every year, for instance looking at one
                  season variability -- will need to fill the options
#    bins -- number of points between 2 chunks
    
    Options:
    MinYears -- minimal year of the time series (default: 1982)
    MaxYears -- Maximal year of the time series (default: 2016)
    str_month, str_day -- month and day the chunk starts 
    end_month, end_day -- month and day the chunk ends

    Output:
    str_id -- index where each chunk starts
    end_id -- index where each chunk ends
    
    '''
    # allocate memory
    str_id = np.zeros(NumChunk)
    end_id = np.zeros(NumChunk)

    if sample==0:
        str_id = range(0,len(dtime)-(tau-1),tau)
        end_id = range(0+tau-1,len(dtime),tau) 

    elif sample==1:
        k=0
        for yy in range(MinYear,MaxYear+1):
            str_chunk = datetime.date(yy,str_month,str_day)
            str_tmp = np.where(dtime == np.datetime64(str_chunk))
            str_id[k] = str_tmp[0][0]
            end_chunk = datetime.date(yy,end_month,end_day)
            end_tmp = np.where(dtime == np.datetime64(end_chunk))
            end_id[k] = end_tmp[0][0]
            k = k + 1

    return str_id, end_id


def PotPred_ZhengFred(ts,tau,NumChunk,str_id,end_id):

    '''
    potential predictability function
    the time series needs to be detrended/deseasoned... 
    before to use the function and depending on the study
    
    The method used here follow the Zheng, Frederiksen, Lou 
    method after:
    Lou et al. 2016 -- doi 10.1007/s00382-016-3229-x
    Frederiksen et al. 2015 -- doi 10.1007/s00382-015-2699-6

    INPUTS:
    ts -- time series to analyse -- 1D
    tau -- length of the chunk
    NumChunk -- total number of chunks
    str_id, end_id -- starting and ending indexes for each chunk

    OUTPUTS
    Var_interC -- inter chunk variance
    Var_noise -- noise component variance
    Var_slow -- slow component variance
    p -- potential predictability ratio following Lou et al. 2016

    Author: Eva Cougnon
    Created: Mar 2017
    Last modified: Apr 2017

    '''

# checking that the time series is not empty
    if ((np.isnan(np.min(ts)) == False) & \
        (np.isnan(np.max(ts)) == False)):
# chunck time series
        ts_chunk=np.empty(tau*NumChunk)
        ts_chunk.fill(np.nan)
# chunk mean allocation
        xc_0=np.empty(NumChunk)
        xc_0.fill(np.nan)
        tmp_A=np.empty(NumChunk) # used for the unpredictable chunk
        tmp_A.fill(np.nan)
        kk=0
        for c in range(0,NumChunk): # run through each chunk
# chunk meam
            xc_0[c] = np.nanmean(ts[int(str_id[c]):int(end_id[c])+1])
            tmp_A[c] = np.nansum((ts[int(str_id[c]+1) \
                                     :int(end_id[c])+1] \
                                  - ts[int(str_id[c]) \
                                       :int(end_id[c])])**2)
            ts_chunk[kk:kk+tau] = ts[int(str_id[c]) \
                                     :int(end_id[c])+1]
            kk = kk + tau
# mean of all the chunks
        xc_0all = np.nanmean(xc_0)
# total inter chunk variance
        tmp=np.zeros(NumChunk)
        for c in range(0,NumChunk):
            tmp[c] = (xc_0[c]-xc_0all)**2
        Var_interC = 1/(NumChunk-1) * np.nansum(tmp)
## estimatingthe variance of the unpredictible chunk 
# estimating the A coeficient
        A = (1/((tau-1)*NumChunk)) * np.nansum(tmp_A)
# remove the chunk length running mean filter of the combined chunk time series
        ts_box=np.empty(len(ts_chunk))
        ts_box.fill(np.nan)
        window = tau
        weights = np.repeat(1.0,window)/window
        ts_box = np.convolve(ts_chunk, weights,'valid')
#        ts_a = ts_chunk[int(tau/2):-int(tau/2)]-ts_box
        ts_a = ts_chunk[int(tau/2):-int(tau/2)+1]-ts_box
# auto correlation coef with lag-1
        coef_corr = np.corrcoef(ts_a[1:],ts_a[0:-1])
        alpha = coef_corr[0][1]
# expectation of the noise (residual between the data and the potentially
# predictable component
        gamma = A/(2-alpha-alpha)
# variance of the noise over all chunks
        tmp=np.zeros(tau)
        for c in range(0,tau-1):
            tmp[c] = (tau-(c+1))*(alpha**(c+1)+alpha**(c+1)) 
        Var_noise = (gamma/(tau**2))*(np.sum(tmp)+tau)
#  variance ofthe potentially predictable component
        Var_slow = Var_interC - Var_noise
# Potential inter-chunk predictability
        p = Var_slow / Var_interC

    else:
        Var_interC = np.nan
        Var_noise = np.nan
        Var_slow = np.nan
        p = np.nan

    return Var_interC, Var_noise, Var_slow, p



def PotPred_vStorchZwiers(ts,tau,NumChunk,str_id,end_id):

    '''
    potential predictability function
    the time series needs to be detrended/deseasoned... 
    before to use the function and depending on the study
    
    The method used here follow the von Storch and Zwiers
    method after their book:
    Statistical Analysis in Climate Research Chapter 17

    INPUTS:
    ts -- time series to analyse
    tau -- length of the chunk
    NumChunk -- total number of chunks
    str_id, end_id -- starting and ending indexes for each chunk

    OUTPUTS
    Var_interC -- inter chunk variance
    Var_noise -- noise component variance
    p -- potential predictability ratio following Lou et al. 2016
    F90 -- significance at 90th percentile assuming the F-distribution
    F95 -- same than above for the 95th

    Author: Eva Cougnon
    Created: Apr 2017
    Last modified: Apr 2017

    '''

# checking that the time series is not empty
    if ((np.isnan(np.min(ts)) == False) & \
        (np.isnan(np.max(ts)) == False)):
# chunk mean allocation
        t_j=np.empty(NumChunk)
        t_j.fill(np.nan)
# periodogram of each chunk
        pxx_j = np.empty([NumChunk,tau])#np.empty([NumChunk,int(tau/2)+1])
        pxx_j.fill(np.nan)# power spectrum of each chunk
        f_j = np.empty([NumChunk,tau])#np.empty([NumChunk,int(tau/2)+1])
        f_j.fill(np.nan)# corresponding frequencies
        for c in range(0,NumChunk): # run through each chunk
# chunk meam
            t_j[c] = np.nanmean(ts[int(str_id[c]):int(end_id[c])+1])
# spectrum of each chunck for the fast varying component
            f_j[c,:],pxx_j[c,:]=signal.periodogram(ts[int(str_id[c]) \
                                                      :int(end_id[c])+1],detrend=False, \
                                                   return_onesided=False, scaling='density')
# mean of all the chunks
        t_jall = np.nanmean(t_j)
# total inter chunk variance
        tmp=np.zeros(NumChunk)
        for c in range(0,NumChunk):
            tmp[c] = (t_j[c]-t_jall)**2
        Var_interC = 1/(NumChunk-1) * np.nansum(tmp)
# inter chunk variance of the fast varying process
# average periodograms in time -- estimate the spectrum for all the chunks
        pxx_jall=np.mean(pxx_j,axis=0)

        idx = np.argmin(np.abs(f_j[0][:] - (1/tau)))
# variance of the fast varying component
        if idx!=0:
            Var_noise = (1/tau) * pxx_jall[idx]*2   # when double sided
        elif idx==0:
            Var_noise = (1/tau) * pxx_jall[idx]

# variance ratio
        p = Var_interC / Var_noise
# testing the null hypothesis to find the significance
        F90 = stats.f.ppf(0.90,tau-1,2*NumChunk)
        F95 = stats.f.ppf(0.95,tau-1,2*NumChunk)

    else:
        Var_interC = np.nan
        Var_noise = np.nan
        p = np.nan
        F90 = np.nan
        F95 = np.nan
 
    return Var_interC, Var_noise, p, F90, F95



def cov_unpred_testing(ts_mtx, tau, NumChunk, str_id, end_id): 
# cov_unpred_testing(ts1, ts2, alpha_1, alpha_2, tau, NumChunk, str_id, end_id):
# cov_unpred_testing(ts_mtx, tau, NumChunk, str_id, end_id):
    ''' 
    estimating the co-variance matrix of the unpredictible component

    alpha calculation should be part of the function
    ts_mtx -- time series with a dimension of (time, location) 
		with land points removed, for instance time can be of 
		length of 12*35 if each chunk is 12 months for 35 years
    tau -- chunk length
    NumChunk -- total number of chunks
    str_id, end_id -- start and end indeces of each chunk

    Author: Eva C.
    Created: May 2017
    Last modif: Jul 2017 -- adding the alpha calculation

    '''

# create the time series of all the chunks combined
    ts_mtx_chunk = np.empty((tau*NumChunk,len(ts_mtx[0,:])))
    ts_mtx_chunk.fill(np.nan)
    kk=0
    for c in range(0,NumChunk): # run through each chunk
        ts_mtx_chunk[kk:kk+tau,:] = ts_mtx[str_id[c]:end_id[c]+1,:]
        kk = kk + tau

# Apply a low pass filter to extract the noise: remove the chunk 
#  length running mean filter of the combined chunk time series 
#  (ts_mtx_chunk)
# NB: moving_average is a function from this file
    ts_box = np.apply_along_axis(moving_average,0,ts_mtx_chunk[:,:],tau)
# remove low pass filter -- keep the noise part  
# depending on tau you may need: int(tau/2):-int(tau/2),:] (works with tau=365)
#	or int(tau/2):-int(tau/2)+1,:] (works with tau =12)
    ts_noise = ts_mtx_chunk[int(tau/2):-int(tau/2),:] - ts_box[:,:]
# auto correlation coef with lag-1 of the noise part (low pass filter
#  removed from the ts_mtx_chunk) -- get a vector for each location
# alpha = np.apply_along_axis(eac.get_autocorr_coef,0,ts_noise[1:,:],ts_noise[:-1,:])
# np.apply_along_axis not working....
# NB: get_autocorr_coef is a function from this file
    alpha=np.empty(len(ts_noise[0,:]))
    alpha.fill(np.nan)
    for i in range(0,len(ts_noise[0,:])): # loop through each location
        alpha[i] = get_autocorr_coef(ts_noise[1:,i],ts_noise[:-1,i])
# covariance matrix from Frederiksen et al 2015 for the unpredictable component
    cov_unp = np.empty((len(alpha), len(alpha)))
    cov_unp.fill(np.nan)
    fact_a = 1/((tau-1)*NumChunk) # factor for 'a' calculation -- equation (4) in 
				  # Frederiksen et al 2015 
    tmp_a = np.empty(NumChunk)
    tmp_a.fill(np.nan)
    tt_b=np.arange(1,tau)
    for i in range(0,len(alpha)):
        for j in range(i,len(alpha)):
    
            beta = (np.nansum( (tau - tt_b)*(alpha[i]**tt_b+alpha[j]**tt_b) ) \
                    + tau) / tau**2
# as in their eq (3) combine with eq (4) -- a        
            for tt in range(0,NumChunk):
                tmp_a[tt] = np.nansum((ts_mtx[str_id[tt]+1:end_id[tt]+1,i] - \
                                       ts_mtx[str_id[tt]:end_id[tt],i]) * \
                                      (ts_mtx[str_id[tt]+1:end_id[tt]+1,j] - \
                                       ts_mtx[str_id[tt]:end_id[tt],j]))
            gamma = (fact_a * np.nansum(tmp_a)) / (2 - alpha[i] - alpha[j])

            cov_unp[i,j] = beta * gamma
            cov_unp[j,i] = beta * gamma

    return cov_unp
















