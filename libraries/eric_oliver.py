'''
   functions from Eric Oliver


'''

import numpy as np
import scipy as sp
from scipy import linalg
from scipy import stats


def trend(x, y, alpha=0.05):
    '''
    Calculates the trend of y given the linear
    independent variable x. Outputs the mean,
    trend, and alpha-level (e.g., 0.05 for 95%)
    confidence limit on the trend.
    returns mean, trend, dtrend_95
   
    Written by Eric Oliver 

    '''
    valid = ~np.isnan(y)
    X = np.array([np.ones(len(x)), x-x.mean()])
    beta = linalg.lstsq(X[:,valid].T, y[valid])[0]
    yhat = np.sum(beta*X.T, axis=1)
    t_stat = stats.t.isf(alpha/2, len(x[valid])-2)
    s = np.sqrt(np.sum((y[valid] - yhat[valid])**2) / (len(x[valid])-2))
    Sxx = np.sum(X[1,valid]**2) - (np.sum(X[1,valid])**2)/len(x[valid]) # np.var(X, axis=1)[1]
    return beta[0], beta[1], t_stat * s / np.sqrt(Sxx)


def deseason_harmonic(dat, K, L):
    '''
    deseasoned_data, season, beta = deseason_harmonic(dat, K, L)

    Subtracts the seasonal cycle (season) from the data (data). Season
    is calculated by fitting K harmonics of the annual cycle (as well as the
    mean) to the data. Assumes on year has L elements (i.e., 365 for daily data,
    73 for pentad data, 52 for weekly data, etc.).
    Outputs the deseasonalized data, the season, and the fitting parameters (beta)

    Handles missing values as np.nan's

    Modification: Eva Cougnon (May 2017) -- return nans if the whole time series
                                            is nans

    Written by Eric Oliver, Dalhousie University, 2007-2011
    Adapted from original MATLAB script on 28 November 2012
    '''

#   Valid (non-NaN) elements
    valid = ~np.isnan(dat)

#   ensure that mat is a matrix and a column vector
    dat = np.mat(dat)
    if dat.shape[1]!=0:
        dat = dat.T

#   length of time series and generate time vector
    n = len(dat)
    time = np.mat(np.arange(1,n+1)/(1.*L))

# checking that the matrix is not entirely nans
    if (valid.any()==True):

#   set up mean and harmonics to fit data
        P = np.mat(np.ones((n,1)))
        for k in range(1,K+1):
            P = np.concatenate((P, np.cos(k*2*np.pi*time.T)), 1)
            P = np.concatenate((P, np.sin(k*2*np.pi*time.T)), 1)

#   Remove seasonal cycle by harmonic regression
        beta = (np.linalg.inv(P[valid,:].T*P[valid,:])*P[valid,:].T)*dat[valid]
        season = P*beta
        dat_ds = dat - season

    elif (valid.all()==False):
        dat_ds = np.empty((np.shape(dat)))
        dat_ds.fill(np.nan)
        season = np.empty((np.shape(dat)))
        season.fill(np.nan)
        beta = np.nan 


    return dat_ds, season, beta


def deseason_harmonic_2D(dat, K, L, detrend=False):
    '''
    deseasoned_data, season, beta = deseason_harmonic_2D(dat, K, L)

    Subtracts the seasonal cycle (season) from the data (data). Season
    is calculated by fitting K harmonics of the annual cycle (as well as the
    mean) to the data. Assumes on year has L elements (i.e., 365 for daily data,
    73 for pentad data, 52 for weekly data, etc.).
    Outputs the deseasonalized data, the season, the trend, and the fitting parameters (beta)

    First dimension of dat must be time dimension.
    Does not handle missing values.

    Written by Eric Oliver, Dalhousie University, 2007-2011
    Adapted from original MATLAB script on 28 November 2012
    '''

#   Valid (non-NaN) elements
    valid = ~np.isnan(dat)

#   ensure that mat is a matrix and a column vector
    dat = np.mat(dat)
    #if dat.shape[1]!=0:
    #    dat = dat.T

#   length of time series and generate time vector
    n = dat.shape[0]
    time = np.mat(np.arange(1,n+1)/(1.*L))

#   set up mean and harmonics to fit data
    P = np.mat(np.ones((n,1)))
    for k in range(1,K+1):
        P = np.concatenate((P, np.cos(k*2*np.pi*time.T)), 1)
        P = np.concatenate((P, np.sin(k*2*np.pi*time.T)), 1)
    if detrend:
        P = np.concatenate((P, time.T-time.mean()), 1)

#   Remove seasonal cycle by harmonic regression
    beta = (np.linalg.inv(P.T*P)*P.T)*dat
    season = P[:,0:K*2+1]*beta[0:K*2+1,:]
    if detrend:
        trend = P[:,-1]*beta[-1,:]
    else:
        trend = 0.*season
    dat_ds = dat - season - trend

    return dat_ds, season, trend, beta


#########################################
# calculating the area of a grid cell.  
# First define the functions latlon2km and dxdy.  
# Then given lat and lon vectors you can get the area 
# of each cell from the dxdy function, 
# and define the area dA as dx*dy.
#
# From Eric Olibve (check email 19 Mar 2018)
########################################

def latlon2km(lon1, lat1, lon2, lat2):
    EARTH_RADIUS = 6378.1
    c = np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) + \
        np.cos(np.radians(lon1-lon2)) * np.cos(np.radians(lat1)) * \
        np.cos(np.radians(lat2))
    d = EARTH_RADIUS * np.arccos(c)
    return d

def dxdy(lon, lat):
    '''
    Takes M+1 length lat and N+1 length lon vectors
    and returns MxN 2D arrays of distances across cells
    in x and y directions
    '''
    X = len(lon)-1
    Y = len(lat)-1
    dx = np.zeros((Y,X))
    dy = np.zeros((Y,X))
    for j in range(dx.shape[0]):
        for i in range(dx.shape[1]):
            dx[j,i] = 1e3 * latlon2km(lon[i+1], lat[j], lon[i], lat[j])
            dy[j,i] = 1e3 * latlon2km(lon[i], lat[j+1], lon[i], lat[j])
    return dx, dy



