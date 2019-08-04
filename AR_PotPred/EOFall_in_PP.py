'''

plot the variance at each location of (EOF_i*PC_i)_pred to see/hae an idea
of  what does the EOF explain of the potential predictability map, in other
words what does each mode (from the EOF analysis of the predictable component)
explain the potential predictability map 

'''
# import libraries
import numpy as np
import xarray as xr

from matplotlib import pyplot as plt
import matplotlib
import mpl_toolkits.basemap as bm

fname_pp = '/home/ecougnon/ana/PotPred/PP_SSTa_monthly_1yr_ZFL_Aus_1deg.nc'
fname_eof = '/home/ecougnon/ana/PotPred/EOF/cov_decomp_1yr_ZFL_EOF_Aus_1deg_TMM_Modes1-8'
#figfile = '/home/ecougnon/ana/PotPred/EOF/EOF1-5var_in_PotPred.png'

# data from the potential predictablity calculation

deg_res = 1 #5
lat_min = -55 #-80 #-20 #-55
lat_max = 10 #80 #20 #10
lat_bin = 4 * deg_res
lon_min = 90 #1 #160 #90
lon_max = 180 #360 #270 #180
lon_bin = 4 * deg_res

lat = xr.open_dataset(fname_pp)['lat']
#lat = lat.sel(lat=slice(lat_min,lat_max,lat_bin))
lon = xr.open_dataset(fname_pp)['lon']
#lon = lon.sel(lon=slice(lon_min,lon_max,lon_bin))
PP = xr.open_dataset(fname_pp)['TMM']
#PP = PP.sel(lat=lat, lon=lon)
var_tot = PP[0,:,:] # variance of the inter chunk (total component)
var_slow = PP[2,:,:] # variance of the slow component
var_pp = PP[-1,:,:] # ratio of the variance of the slow over the total component

# data from EOF analysis
data = np.load(fname_eof + '.npz')
MODES_MAP = data['MODES_MAP'].item()
EIGVAL = data['EIGVAL'].item() # to load a dictionary
PC = data['PC'].item()
lat_eof = data['lat']
lon_eof = data['lon']
time_ts = data['time_ts']

NumEOF =4 

PC_id = ['PC1_pred','PC2_pred','PC3_pred','PC4_pred'] #,'PC5_pred']#, \
#         'PC6_pred','PC7_pred','PC8_pred']

var_eof_all = np.empty((NumEOF,np.size(PP,1),np.size(PP,2)))
var_eof_all.fill(np.nan)
for e in range(0,NumEOF):
    var_eof = MODES_MAP['pred_modes_map'][-(e+1),:,:]
    var_PC = PC[PC_id[e]].copy()
    mode_map = np.zeros((len(time_ts),len(var_eof[:,0]),len(var_eof[0,:])))
    for tt in np.arange(0,len(time_ts)):
        mode_map[tt,:,:] = var_eof[:,:]*var_PC[tt]
# calc the normal variance
    mode_map_var = np.var(mode_map,axis=0)
# calc following Frederiksen et al 2016 for the varaince calc of ther
# total component using chunks
# tau -- number of months per chunk
    tau = 12 # 12 months
# time period
    MinYear = 1982
    MaxYear = 2016
    NumYears = MaxYear - MinYear+1
    NumMonths = NumYears*12
    Var_interC = np.empty((len(lat_eof),len(lon_eof)))
    Var_interC.fill(np.nan)
# define the start/end indexes for each chunk
    str_id = range(0,NumMonths,tau)
    end_id = range(tau-1,NumMonths+1,tau)
    NumChunk = len(str_id)
    icnt_id = 0
    for i in range(0,len(lon_eof)):
        jcnt_id = 0
        for j in range(0,len(lat_eof)):
            ts = mode_map[:,j,i]
# checking that the time series is not empty
            if ((np.isnan(np.min(ts)) == False) & \
                (np.isnan(np.max(ts)) == False)):
# chunck time series
                ts_chunk=np.empty(tau*NumChunk)
                ts_chunk.fill(np.nan)
# chunkmean allocation
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
                Var_interC[jcnt_id,icnt_id] = 1/(NumChunk-1) * np.nansum(tmp)
            jcnt_id = jcnt_id + 1
        icnt_id = icnt_id + 1

# ratio that we want to plot (EOF variance over variance of the signal from
# the potential predictability calculation
#var_ratio = mode_map_var/var_tot
    var_eof_all[e,:,:] = Var_interC #/var_tot

var_ratio_sum = np.nansum(var_eof_all,axis=0)/var_tot

## plotting 
domain = [int(lat_min), int(lon_min), int(lat_max), int(lon_max)]
domain_draw = [int(lat_min), int(lon_min), int(lat_max), int(lon_max)]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'
lev_pp = np.arange(0,1+0.1,0.1)
#lev_diff = np.hstack((np.arange(-1,-0.1+0.1,0.1), np.arange(0.1,1+0.1,0.1)))
lev = np.arange(0,100+1,10)
lev_ = np.arange(0,100+1,10)

# PotPred ratio + var(EOF mode)/var(tot) used in PotPred to see how much the
# EOF modes oft he predictable component explain the PotPred ratio 
plt.figure(figsize=(7,16)) 
plt.clf()
ax1 = plt.subplot(3,1,1)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, var_pp, levels= lev_pp, cmap=plt.cm.afmhot_r)
cb=plt.colorbar(ticks = lev_pp, shrink=0.9)
plt.title('Potential predictability ratio on TMM', fontsize=16, y=1.04)

ax2 = plt.subplot(3,1,2)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, var_ratio_sum, levels= lev_pp, cmap=plt.cm.afmhot_r)
cb=plt.colorbar(ticks = lev_pp, shrink=0.9)
plt.title('VAR(mode1-5_pred)/VAR(tot)', fontsize=16, y=1.04)


ax3 = plt.subplot(3,1,3)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, (var_ratio_sum/var_pp)*100, levels= lev, \
             cmap=plt.cm.afmhot_r) #seismic)
# (np.nansum(var_eof_all, axis=0)/var_slow)*100
# (var_ratio_sum/var_pp)*100
#  (1-(var_pp-var_ratio_sum))*100
cb=plt.colorbar(ticks = lev_, shrink=0.9)
plt.title('% of PP explained by the 5 leading modes', fontsize=16, y=1.04)

#plt.savefig(figfile, bbox_inches='tight', dpi=300)

plt.show()










