'''
    plot maps of the EOF analysis calculated with the covariance
    decomposition function (cov_decomp_testing.py or cov_decomp.py)
	also plot the reconstructed time series of each modes for each 
	component and the eigenvalue spectrum to cheack that each
	mode is independent from the other

  Author: Eva C.
  Created: July 2017
  Last Modif:

'''

# import libraries
import numpy as np
import xarray as xr
from scipy.stats import pearsonr

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm


# Load data
fname = '../../ana/PotPred/EOF/cov_decomp_1yr_ZFL_SVD_Aus_daily_1deg_Modes1-8'
#cov_decomp_1yr_ZFL_SVD_Aus_TMM_trend_1deg_Modes1-8'
#cov_decomp_1yr_ZFL_SVD_Aus_1deg_TMM_Modes1-8'
fname_ppr = '../../ana/PotPred/PP_SSTa_monthly_1yr_ZFL_Aus.nc'
#cov_decomp_1yr_ZFL_SVD_Aus_TMM_trend_1deg_Modes1-8'
#fname_tmp = '/home/ecougnon/ana/PotPred/EOF/cov_decomp_1yr_ZFL_EOF_global_5deg_Tp90'
figfile ='/home/ecougnon/Desktop/WorkingFigures/EOF5-8_tot_Aus_TMM_daily_svd.png'
#figfile2 ='/home/ecougnon/ana/PotPred/EOF/WithSeason/Eigenvalue_spectrum_Aus_TMM_trend.png'
data = np.load(fname + '.npz')
MODES_MAP = data['MODES_MAP'].item()
EIGVAL = data['EIGVAL'].item() # to load a dictionary
PC = data['PC'].item()
data_tmp = np.load(fname + '.npz')
lat = data_tmp['lat']
lon = data_tmp['lon']
time_ts = data_tmp['tim'] #e_ts']
MinYear = 1982
MaxYear = 2016
NumYears = MaxYear-MinYear+1

# what to plot
ID_mode = 'tot_modes_map'
ID_var = 'VAR_tot'
idx_mode = 4
var_eof1_map = MODES_MAP[ID_mode][idx_mode,:,:]
var_PC = PC['PC1_pred'].copy()
var_eof2_map = MODES_MAP[ID_mode][idx_mode+1,:,:]
var_eof3_map = MODES_MAP[ID_mode][idx_mode+2,:,:]
var_eof4_map = MODES_MAP[ID_mode][idx_mode+3,:,:]
var_exp1 = EIGVAL[ID_var][idx_mode]*100 #eigvalP_unp
var_exp2 = EIGVAL[ID_var][idx_mode+1]*100
var_exp3 = EIGVAL[ID_var][idx_mode+2]*100
var_exp4 = EIGVAL[ID_var][idx_mode+3]*100

# load ppr
lat_map = xr.open_dataset(fname_ppr)['lat']
lon_map = xr.open_dataset(fname_ppr)['lon']
PP = xr.open_dataset(fname_ppr)['TMM'] 

'''
## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 
file_enso = np.genfromtxt('/home/data/index/enso/nino34_anomaly.txt', \
			  skip_header=1, skip_footer = 7) 
#skip_footer = 7) # for nino34_anomaly
#MEI                          skip_header=10, skip_footer = 30, delimiter='\t')
#file_other = np.genfromtxt('/home/ecougnon/data/atlantic/amm_sst.txt', \
#                           skip_header=1, skip_footer = 6)

str_id = np.nonzero((file_enso[:,0]>(MinYear-1)) \
                    & (file_enso[:,0]<(MinYear+1)))
enso_monthly_tmp = np.empty(NumYears*12)
#other_monthly = np.empty(NumYears*12)
k=0
for yy in np.arange(str_id[0][0],len(file_enso[:,0])):
    for mm in np.arange(1,12+1):
        enso_monthly_tmp[k] = file_enso[yy,mm]
#        other_monthly[k] = file_other[yy,mm]
        k = k + 1
# standardise the index
enso_monthly = (enso_monthly_tmp-np.nanmean(enso_monthly_tmp)) \
               /np.nanstd(enso_monthly_tmp)
'''

#######################################
lat_min = np.min(lat)
lat_max = np.max(lat)
lon_min = np.min(lon)
lon_max = np.max(lon)
'''
lat_bin = 4 * 5
lon_bin = 4 * 5
fname2 = '/home/ecougnon/ana/SSTa_monthly_extremes_global.nc'
lat = xr.open_dataset(fname2)['lat']
lat = lat.sel(lat=slice(lat_min,lat_max,lat_bin))
lon = xr.open_dataset(fname2)['lon']
lon = lon.sel(lon=slice(lon_min,lon_max,lon_bin))
time_ts = xr.open_dataset(fname2)['time']
'''
############################################

## plotting 
'''
domain = [int(lat_min)-1, int(lon_min), int(lat_max)+1, int(lon_max)+1] 
domain_draw = [int(lat_min)-1, int(lon_min), int(lat_max)+1, int(lon_max)+1] 
'''
domain = [-55, 90, 10, 180] #[-80, -180, 85, 180] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon, lat)
llon_, llat_ = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'
lev = np.hstack((np.arange(-0.07,-0.01+0.01,0.02), np.arange(0.01,0.07+0.01,0.02)))

# plotting the 3 first EOF on the same figure

plt.figure(figsize=(17,4)) #5,12))
plt.clf()
ax1 = plt.subplot(1,4,1) #3,1,1)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
lonproj_, latproj_ = proj(llon_, llat_)
plt.contourf(lonproj, latproj, var_eof1_map, levels= lev, cmap=plt.cm.seismic) 
cb=plt.colorbar(ticks=lev,shrink=0.7)
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.75], colors='C2')
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.2], colors='C1')
plt.title('EOF' + str(idx_mode+1) +', ' + str(round(var_exp1,1)) +'%', fontsize=16, y=1.04)

ax2 = plt.subplot(1,4,2) #3,1,2)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
lonproj_, latproj_ = proj(llon_, llat_)
plt.contourf(lonproj, latproj, var_eof2_map, levels= lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.75], colors='C2')
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.2], colors='C1')
plt.title('EOF' + str(idx_mode+2) +', ' + str(round(var_exp2,1)) +'%', fontsize=16, y=1.04)

ax3 = plt.subplot(1,4,3) #3,1,3)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
lonproj_, latproj_ = proj(llon_, llat_)
plt.contourf(lonproj, latproj, var_eof3_map, levels= lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.75], colors='C2')
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.2], colors='C1')
plt.title('EOF' + str(idx_mode+3) +', ' + str(round(var_exp3,1)) +'%', fontsize=16, y=1.04)

ax4 = plt.subplot(1,4,4)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
lonproj_, latproj_ = proj(llon_, llat_)
plt.contourf(lonproj, latproj, var_eof4_map, levels= lev, cmap=plt.cm.seismic)
cb=plt.colorbar(ticks=lev,shrink=0.7)
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.75], colors='C2')
plt.contour(lonproj_, latproj_, PP[3,:,:], levels=[0.2], colors='C1')
plt.title('EOF' + str(idx_mode+4) +', ' + str(round(var_exp4,1)) +'%', fontsize=16, y=1.04)


plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)
#plt.show(block=False)


'''
# EOF map with their corresponding PCs

plt.figure(figsize=(7,11)) #12,11)) #7,11)) 
plt.clf()
ax1 = plt.subplot(2,1,1)
#proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, var_eof_map, levels= lev, cmap=plt.cm.seismic) 
cb=plt.colorbar(ticks=lev,shrink=0.9)
#cb.ax1.tick_params(labelsize=14)
plt.title('EOF 1, % -- pred (1x1) OTEp10', fontsize=16, y=1.04)

ax2 = plt.subplot(2,1,2)
## plot the reconstructed time series for each mode
var_comp = enso_monthly*10 # factor 10 to be readable on the figure
r1, _ = pearsonr(var_PC,var_comp)
ax_pc = plt.gca()
plt.plot(time_ts,var_PC)
plt.plot(time_ts,var_comp)
plt.xlim(np.min(time_ts), np.max(time_ts))
ax_pc.annotate('r = {:.2f}'.format(r1), xy=(.1,.95), xycoords=ax_pc.transAxes)
plt.title('Corresponding PC with NINO3.4 index')
plt.grid()
plt.savefig(figfile, bbox_inches='tight', dpi=300)
#plt.show()
'''


'''
## eigenvalues spectrum -- North et al. 1982
plt.figure()
plt.errorbar(np.arange(0,8,1),EIGVAL['VAR_tot'][:8], \
             yerr=EIGVAL['delta_VAR_tot'][:8],fmt='.')
plt.errorbar(np.arange(0,8,1),EIGVAL['VAR_pred'][:8], \
             yerr=EIGVAL['delta_VAR_pred'][:8],fmt='.')
plt.errorbar(np.arange(0,8,1),EIGVAL['VAR_unp'][:8], \
             yerr=EIGVAL['delta_VAR_unp'][:8],fmt='.')
plt.legend(['tot','pred','unp'])
plt.title('Eigenvalues spectrum -- North et al. 1982')
plt.xlabel('EOF modes')
plt.ylabel('Eigenvalue (%)')
plt.grid()
plt.savefig(figfile2, bbox_inches='tight', dpi=300)
plt.show()
'''





