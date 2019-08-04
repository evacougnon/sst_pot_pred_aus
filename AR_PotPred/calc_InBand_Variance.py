'''
        Use signal.periodogram to estimate an appropriate
        time scale to apply with the PotPred methods
 
        Created: Aug 2018
        Author: Eva C.
        Last Modif.: Oct 2018 (moved the danbpower function to eac_useful.py
'''

##################
# import libraries
###################
import numpy as np
import xarray as xr
from scipy import signal
#import dask.array as da
#import bottleneck as bn

#from xarray.core import dask_array_ops

import sys
sys.path.insert(0,'../libraries/')
import eac_useful as eac
import eric_oliver as eo

import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import cmocean

###########################
# I/O file names
#########################
#outfile = 
#figfile = '/home/ecougnon/ana/PotPred/InBand/'

############################
# define region and load data
# specify a lat and lon bin if subsetting!
# TMM: monthly temperature 
# Tp90: monthly 90th percentile
# TMX: monthly max
# Tp10: monthly 10th percentile
# TMN: monthly min
#
# 'HadISST' to use for the HadISST file
# already monthly SSTs
############################
WHICH = 'TMM' #'HadISST' #'TMM'
lat_min = -56
lat_max = 10
lat_bin = 1
lon_min = 90
lon_max = 180
lon_bin = 1
lat_wanted = -5
lon_wanted = 170
if WHICH == 'TMM':
    fname = '../tmp_data2plot/SSTa_monthly_extremes_Aus.nc'
#?SSTa_daily_Aus.nc'
#SSTa_monthly_extremes_Aus.nc'
    lat = xr.open_dataset(fname)['lat'].sel(lat=slice(lat_min,lat_max,lat_bin))
    lon = xr.open_dataset(fname)['lon'].sel(lon=slice(lon_min,lon_max,lon_bin))
    tim = xr.open_dataset(fname)['time']
    SSTa_TMm = xr.open_dataset(fname)[WHICH].sel(time=tim, lat=lat, lon=lon)
    X = len(lon)
    Y = len(lat)
elif WHICH == 'HadISST':
    fname = '/home/ecougnon/Desktop/WorkingFigures/HadSST3/HadISST_sst.nc'
    lat = xr.open_dataset(fname)['latitude']. \
             sel(latitude=slice(lat_max,lat_min,lat_bin))
    lon = xr.open_dataset(fname)['longitude']. \
             sel(longitude=slice(lon_min,lon_max,lon_bin))
    tim = xr.open_dataset(fname)['time'].sel(time=slice('1960-01-01','2016-12-31'))
    SST = xr.open_dataset(fname)['sst'].sel(time=tim, latitude=lat, longitude=lon)
    X = len(lon)
    Y = len(lat)
# detrend the time series
    SST=np.array(SST)
    SST[SST<-100]=np.nan
    SST_flat = SST.reshape(len(tim),X*Y)
    dsst = np.empty((len(tim),X*Y))
    dsst.fill(np.nan)
    tt=np.arange(0,len(tim))
    for i in range(0,len(SST_flat[0,:])):
        valid = ~np.isnan(SST_flat[:,i])
        if (valid.any()==True):
            y = SST_flat[:,i]
            mean, trend, alpha = eo.trend(tt,y)
            dsst[:,i] = y - (tt*trend) -mean
#            dsst[:,i] = signal.detrend(SST_flat[valid,i], axis=0, type='linear')
        elif (valid.all()==False):
            dsst[:,i] = np.nan
        tmp_dsea, sea, beta = eo.deseason_harmonic(dsst[:,i],4,12)
        dsst[:,i] = np.squeeze(tmp_dsea[:,0])
    SSTa_TMm = dsst.reshape(len(tim),Y,X)

#figfile ='/v_Munk_Drive/ecougnon/ana/InBand_Variance/InBandVar_SmoothPeriodogram_Scaling2TotVar_1-3-7-10.png'
#figfile_ ='/v_Munk_Drive/ecougnon/ana/InBand_Variance/InBandVar_SmoothPeriodogram_Scaling2TotVar_1-3-7-10.eps'

figfile_ = '../tmp_data2plot/InBandVar_SmoothPeriodogram_Scaling2TotVar_1-3-7-10.eps'

'''
# testing at 1 location
#test = np.squeeze(SSTa_TMm[:,0,0])
test =  xr.open_dataset(fname)[WHICH].sel(lat=lat_wanted, lon=lon_wanted, \
                                          method='nearest')
#sampling rate -- mmonthy
fs = 1#2 #/((365/12)*24*3600) # cycle per day (1/(365/12))
# second per month: (365/12)*24*3600
f_, Pxx_ = signal.periodogram(test,fs, window='hanning', \
                              detrend=False, scaling='density')
# f_1(below) is the same than a normal periodogram 
#f_1, Pxx_den_1 = signal.welch(test, fs, window='hanning', nperseg=len(test), \
#                              detrend=False, scaling='density')
# trying different size of segments and overlap len(test)/2 = 2 segments
f_2, Pxx_den_2 = signal.welch(test, fs, window='hanning', nperseg=len(test)/2, \
                              noverlap=0, detrend=False, scaling='density')
#f_4no, Pxx_den_4no = signal.welch(test, fs, window='hanning', nperseg=len(test)/35, \
#                         noverlap=0, detrend=False, scaling='density')

#################
# plotting one periodogram
#################
plt.figure()
plt.plot(f_,Pxx_)
plt.plot(f_2,Pxx_den_2)
plt.grid()
plt.xlabel('cycle/month')
plt.ylabel('Power spectrum density')
plt.legend(['periodogram','welch number of segm=2'])
plt.title('NW -- 20:30S - 105:115E')
plt.xlim([0, 0.2])
plt.show()
'''

####################################
# calc the power spectrum density
###################################

var_tot = np.var(SSTa_TMm,axis=0)
min_tot = 0 #np.nanmin(var_tot)
max_tot = np.nanmax(var_tot)
#sampling rate -- monthy
fs = 1 #/((365/12)*24*3600) # cycle per day (1/(365/12))
# second per month: (365/12)*24*3600

f_per, pxx_per = np.apply_along_axis(signal.periodogram,0,SSTa_TMm,fs, \
                                     window ='hanning',detrend=False, \
                                     scaling='density')
f_pwe, pxx_pwe = np.apply_along_axis(signal.welch,0,SSTa_TMm,fs, \
                                     window ='hanning',nperseg=len(tim)/2, \
                                     noverlap=0, detrend=False, \
                                     scaling='density')
f_dblside, Pxx_dblside = np.apply_along_axis(signal.periodogram,0,SSTa_TMm,fs, \
                                         window ='hanning',detrend=False, \
                                         return_onesided=False, scaling='density')
# ressorting index for the double sided
f_dblside_ = np.nan*np.zeros((len(f_dblside[:,0,0]),Y,X))
Pxx_dblside_ = np.nan*np.zeros((len(f_dblside[:,0,0]),Y,X))
Pxx_dblside_smooth = np.nan*np.zeros((len(f_dblside[:,0,0]),Y,X))
window = 3
for ii in range(0,Y):
    for jj in range(0,X):
        id_sort = np.argsort(f_dblside[:,ii,jj],axis=0)
        f_dblside_[:,ii,jj] = f_dblside[id_sort,ii,jj]
        Pxx_dblside_[:,ii,jj] = Pxx_dblside[id_sort,ii,jj]*2
        Pxx_dblside_[int(len(Pxx_dblside_)/2),ii,jj] = \
                    Pxx_dblside_[int(len(Pxx_dblside_)/2),ii,jj]/2
        Pxx_dblside_smooth[int(window/2):-int(window/2),ii,jj]=  \
                           eac.movingaverage(Pxx_dblside_[:,ii,jj],window)


InBand0_p = np.nan*np.zeros((Y,X)) # band between <1 years
InBand1_p = np.nan*np.zeros((Y,X)) # band between 1-2 years or 1-3
#InBand2_p = np.nan*np.zeros((Y,X)) # band between 2-4 years
InBand4_p = np.nan*np.zeros((Y,X)) # band between 4-7 years or 3-7
InBand7_p = np.nan*np.zeros((Y,X)) # band between 7-10 years
InBand10_p = np.nan*np.zeros((Y,X)) # band between >10 years
InBand_all_p = np.nan*np.zeros((Y,X)) # band between all
pxx = Pxx_dblside_smooth[int(window/2):-int(window/2),:,:] #pxx_pwe #pxx_per
f = f_dblside_[int(window/2):-int(window/2),:,:] #f_pwe #f_per
max_f = np.max(np.max(f,axis=0))
for ii in range(0,Y):
    for jj in range(0,X):
        InBand0_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0.08,max_f)
        InBand1_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0.0278,0.08)
#0.042,0.08)
#        InBand2_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0.021,0.042)
        InBand4_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0.012,0.0278)
#0.012,0.021)
        InBand7_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0.0083,0.012)
        InBand10_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0,0.0083)
        InBand_all_p[ii,jj] = eac.bandpower_Pxx(pxx[:,ii,jj],f[:,ii,jj],0,max_f)

#################################
# Scale the data with min max of 
# the in-band variance
#################################
min_0 = np.nanmin(InBand0_p)
max_0 = np.nanmax(InBand0_p)
min_1 = np.nanmin(InBand1_p)
max_1 = np.nanmax(InBand1_p)
#min_2 = np.nanmin(InBand2_p)
#max_2 = np.nanmax(InBand2_p)
min_4 = np.nanmin(InBand4_p)
max_4 = np.nanmax(InBand4_p)
min_7 = np.nanmin(InBand7_p)
max_7 = np.nanmax(InBand7_p)
min_10 = np.nanmin(InBand10_p)
max_10 = np.nanmax(InBand10_p)
min_all = 0 #np.nanmin(InBand_all_p)
max_all = np.nanmax(InBand_all_p)
# normalised by the total variance
InBand0_p_norm = 100 * (InBand0_p/InBand_all_p)
#((InBand0_p - min_all)/(max_all-min_all))
InBand1_p_norm = 100 * (InBand1_p/InBand_all_p)
#((InBand1_p - min_all)/ \(max_all-min_all))
#InBand2_p_norm = 100 * ((InBand2_p - min_all)/ \
#                       (max_all-min_all))
InBand4_p_norm = 100 * (InBand4_p/InBand_all_p)
#((InBand4_p - min_all)/ \(max_all-min_all))
InBand7_p_norm = 100 * (InBand7_p/InBand_all_p)
#((InBand7_p - min_all)/ \(max_all-min_all))
InBand10_p_norm = 100 * (InBand10_p/InBand_all_p)
#((InBand10_p - min_all)/ \(max_all-min_all))
# normalised by its own min/max
InBand0_p_norm_ = 100 * ((InBand0_p - min_0)/ \
                       (max_0-min_0))
InBand1_p_norm_ = 100 * ((InBand1_p - min_1)/ \
                       (max_1-min_1))
#InBand2_p_norm_ = 100 * ((InBand2_p - min_2)/ \
#                       (max_2-min_2))
InBand4_p_norm_ = 100 * ((InBand4_p - min_4)/ \
                       (max_4-min_4))
InBand7_p_norm_ = 100 * ((InBand7_p - min_7)/ \
                       (max_7-min_7))
InBand10_p_norm_ = 100 * ((InBand10_p - min_10)/ \
                       (max_10-min_10))


#################################
# plotting
#################################
levs=np.arange(0,100+10,10)

plt.figure(figsize=(17,8))
ax1=plt.subplot(2,3,3,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#plt.contourf(InBand1_p, cmap=plt.cm.viridis)
#plt.contourf(InBand1_p_norm_,levels=levs, cmap=plt.cm.viridis)
#plt.contourf(lon, lat, InBand1_p_norm,levels=np.arange(np.nanmin(InBand1_p_norm), \
#             np.nanmax(InBand1_p_norm)), cmap=plt.cm.afmhot_r)
plt.contourf(lon, lat, InBand1_p_norm,levels=np.arange(5,55+5,5),cmap=plt.cm.afmhot_r)
levels=np.arange(-2,2.2,0.2),
#plt.gca().invert_yaxis() # to use with HAdISST dataset
cb=plt.colorbar(ticks=np.arange(5,55+5,10),shrink=0.7)
cb.ax.tick_params(labelsize=10)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('1-3yrs; ' + str(max(0,round(np.nanmin(InBand1_p_norm),1))) +' - ' + \
           str(round(np.nanmax(InBand1_p_norm),1)) +' %')
'''
plt.subplot(2,3,3)
plt.contourf(InBand2_p, cmap=plt.cm.viridis)
#plt.contourf(InBand2_p_norm_,levels=levs, cmap=plt.cm.viridis)
#plt.contourf(InBand2_p_norm,levels=np.arange(np.nanmin(InBand2_p_norm),np.nanmax(InBand2_p_norm)), cmap=plt.cm.viridis)
plt.gca().invert_yaxis()
plt.colorbar()
plt.title('2-4yrs; ' + str(max(0,round(np.nanmin(InBand2_p_norm),1))) +' - ' + \
           str(round(np.nanmax(InBand2_p_norm),1)) +' %')
'''
ax1=plt.subplot(2,3,4,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#plt.contourf(InBand4_p, cmap=plt.cm.viridis)
#plt.contourf(InBand4_p_norm_,levels=levs, cmap=plt.cm.viridis)
#plt.contourf(lon, lat, InBand4_p_norm,levels=np.arange(np.nanmin(InBand4_p_norm), \
#             np.nanmax(InBand4_p_norm)), cmap=plt.cm.afmhot_r)
plt.contourf(lon, lat, InBand4_p_norm,levels=np.arange(0,30+5,5), cmap=plt.cm.afmhot_r)
#plt.gca().invert_yaxis()
cb=plt.colorbar(ticks=np.arange(0,30+5,5),shrink=0.7)
cb.ax.tick_params(labelsize=10)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('3-7yrs; ' + str(max(0,round(np.nanmin(InBand4_p_norm),1))) +' - ' + \
           str(round(np.nanmax(InBand4_p_norm),1)) +' %')

ax1=plt.subplot(2,3,5,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#plt.contourf(InBand7_p, cmap=plt.cm.viridis)
#plt.contourf(InBand7_p_norm_,levels=levs, cmap=plt.cm.viridis)
#plt.contourf(lon, lat, InBand7_p_norm,levels=np.arange(np.nanmin(InBand7_p_norm), \
#             np.nanmax(InBand7_p_norm)), cmap=plt.cm.afmhot_r)
plt.contourf(lon, lat, InBand7_p_norm,levels=np.arange(0,30+5,5), cmap=plt.cm.afmhot_r)
#plt.gca().invert_yaxis()
cb=plt.colorbar(ticks=np.arange(0,30+5,5),shrink=0.7)
cb.ax.tick_params(labelsize=10)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('7-10yrs; ' + str(max(0,round(np.nanmin(InBand7_p_norm),1))) +' - ' + \
           str(round(np.nanmax(InBand7_p_norm),1)) +' %')

ax1=plt.subplot(2,3,2,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#plt.contourf(InBand0_p, cmap=plt.cm.viridis)
#plt.contourf(InBand0_p_norm_,levels=levs, cmap=plt.cm.viridis)
#plt.contourf(lon, lat, InBand0_p_norm,levels=np.arange(np.nanmin(InBand0_p_norm), \
#             np.nanmax(InBand0_p_norm)), cmap=plt.cm.afmhot_r)
plt.contourf(lon, lat, InBand0_p_norm,levels=np.arange(10,80+5,5), cmap=plt.cm.afmhot_r)
#plt.gca().invert_yaxis()
cb=plt.colorbar(ticks=np.arange(10,80+5,10),shrink=0.7)
cb.ax.tick_params(labelsize=10)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('<1yrs; ' + str(max(0,round(np.nanmin(InBand0_p_norm),1))) +' - ' + \
           str(round(np.nanmax(InBand0_p_norm),1)) +' %')

ax1=plt.subplot(2,3,6,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#plt.contourf(InBand10_p, cmap=plt.cm.viridis)
#plt.contourf(InBand10_p_norm_,levels=levs, cmap=plt.cm.viridis)
#plt.contourf(lon, lat, InBand10_p_norm,levels=np.arange(np.nanmin(InBand10_p_norm), \
#             np.nanmax(InBand10_p_norm)), cmap=plt.cm.afmhot_r)
plt.contourf(lon, lat, InBand10_p_norm,levels=np.arange(0,30+5,5), cmap=plt.cm.afmhot_r)
#plt.gca().invert_yaxis()
cb=plt.colorbar(ticks=np.arange(0,30+5,5),shrink=0.7)
cb.ax.tick_params(labelsize=10)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('>10yrs; ' + str(max(0,round(np.nanmin(InBand10_p_norm),1))) +' - ' + \
           str(round(np.nanmax(InBand10_p_norm),1)) +' %')

ax1=plt.subplot(2,3,1,projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND)
ax1.coastlines('50m', linewidth=0.8)
ax1.gridlines()
#plt.contourf(InBand10_p, cmap=plt.cm.viridis)
#plt.contourf(InBand10_p_norm_,levels=levs, cmap=plt.cm.viridis)
plt.contourf(lon, lat, InBand_all_p,levels=np.arange(0, 1.2, 0.1), \
             cmap=plt.cm.afmhot_r)
#plt.gca().invert_yaxis()
cb=plt.colorbar(shrink=0.7)
cb.ax.tick_params(labelsize=10)
ax1.set_xlim([90, 180])
ax1.set_ylim([-55, 10])
ax1.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
plt.title('Total Variance')

#plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)
plt.savefig(figfile_, bbox_inches='tight', format='eps', dpi=300)

'''
plt.figure()
#plt.plot(f_per[:,92,97],pxx_per[:,92,97])
#plt.plot(f_pwe[:,92,97],pxx_pwe[:,92,97])
#plt.plot(f_per[:,29,239],pxx_per[:,29,239])
#plt.plot(f_pwe[:,29,239],pxx_pwe[:,29,239])
#plt.legend(['periodogram','pwelch'])
#plt.plot(f_dblside_[:,29,239],Pxx_dblside_[:,29,239],'k-*')
#plt.plot(f_dblside_[:,29,239],Pxx_dblside_smooth[:,29,239],'-*')
plt.plot(f_dblside_[:,55,57],Pxx_dblside_[:,55,57],'k-*')
plt.plot(f_dblside_[:,55,57],Pxx_dblside_smooth[:,55,57],'-*')
plt.legend(['not smoothed','smoothed'])
plt.plot([0.012, 0.012],[0, 100],'k--') #7 yr
plt.plot([0.0167, 0.0167],[0, 100],'k--') #5yr
#plt.plot([0.021, 0.021],[0, 100],'k--') # 4yr
plt.plot([0.0278, 0.0278],[0, 100],'k--') # 3yr
#plt.plot([0.042, 0.042],[0, 100],'k--') # 2yr
plt.plot([0.08, 0.08],[0, 100],'k--') # 1yr
plt.plot([0.0083, 0.0083],[0, 100],'k--',alpha=0.3)
plt.plot([0.167, 0.167],[0, 100],'k--',alpha=0.3)
plt.ylim([0,14])
plt.xlim([-0.5,0.5])
plt.title('lat: -48, lon: 150; with 10, 7, 4, 2, 1yr and 6 months marks')
plt.grid()
'''

plt.show(block=False)

