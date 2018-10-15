'''

    read and plot outputs from the potential predictability 
    analysis

'''


# import libraries
import numpy as np
import xarray as xr

#from matplotlib import pyplot as plt
#import matplotlib
#import mpl_toolkits.basemap as bm

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


fname = '../../ana/PotPred/PP_SSTa_daily_1yr_vSZ_Aus.nc'
fname_='../../ana/PotPred/PP_SSTa_monthly_1yr_vSZ_Aus.nc'
#'../../ana/PotPred/PP_SSTa_monthly_1yr_vSZ_Aus.nc'
#fname = '../../ana/PotPred/PP_SSTa_monthly_1yr_ZFL_Aus.nc'
#PP_SSTa_trend_monthly_1yr_ZFL_Aus.nc'
#PP_35yrsStart_monthly_1yr_ZFL_Aus.nc'
#PP_SSTa_monthly_1yr_ZFL_Aus.nc'
#PP_OTEs_NumDaysP90_monthly_1yr_ZFL_Aus.nc'
#PP_SSTa_monthly_with_season_1yr_ZFL_Aus.nc'
#fname = '/home/ecougnon/ana/PotPred/PP_HadISST_monthly_5yr_ZFL_global.nc'
#figfile = '/home/ecougnon/ana/PotPred/PPratio_SSTa_TMM_1yr_ZFL_sign_vSZ.png'
#PPratio_OTE_p90_1yr_ZFL_Aus.png'
figfile = '/home/ecougnon/Desktop/WorkingFigures/vSZ/PPR_vSZ_day.png'

# when using a nc file
lat_map = xr.open_dataset(fname)['lat']
lon_map = xr.open_dataset(fname)['lon']
PP = xr.open_dataset(fname)['TMM'] #Pdays_p90'] #TMM']
PP=PP.transpose('PP_keys','lat','lon')
PP_ = xr.open_dataset(fname_)['TMM']
'''
data = np.load(fname+'.npz')
lon_map = data['lon_map']
lat_map = data['lat_map']
PP1 = data['PP_Tp90'].item()
PP2 = data['PP_TMM'].item()
'''
# plot setting
'''
domain = [-55, 90, 10, 180] #[-80, 0, 85, 360] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-80, 0, 85, 360] #[-50, 90, 10, 180]
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'
'''

ax=plt.figure(figsize=(10,8)) #12,6)) # (10,8)
plt.clf()
#proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
#ax.add_feature(cfeature.COASTLINE)
ax.coastlines('50m', linewidth=0.8)
ax.gridlines()
#proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
#                    zorder=None, alpha=None)
#proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
#                   labels=[True,False,False,False], fontsize=14)
#proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
#                   labels=[False,False,False,True], fontsize=14)
#lonproj, latproj = proj(llon, llat)
plt.contourf(lon_map, lat_map, ((PP[0,:,:]- PP[1,:,:]) / PP[0,:,:]), \
             levels=np.arange(0,1.1,0.1), cmap=plt.cm.afmhot_r, \
             transform=ccrs.PlateCarree())
#plt.contourf(lonproj, latproj, ((PP[0,:,:]- PP[1,:,:]) / PP[0,:,:]), \
#levels=np.arange(0,3+0.1,0.1), \
#                  cmap=plt.cm.afmhot_r)  #Oranges)#YlOrBr)
#cb=plt.colorbar(ticks=np.arange(0,1+0.1,0.1),shrink=0.9)
#plt.contourf(lonproj, latproj, PP1['p'] - PP2['p'], levels=np.arange(-0.5,0.5+0.01,0.01), \
#             cmap=plt.cm.seismic) #BrBG)
cb=plt.colorbar() #ticks=np.arange(-0.5,0.5+0.01,0.1),shrink=0.9)
cb.ax.tick_params(labelsize=14)
ax.set_xlim([90, 180])
ax.set_ylim([-55, 10])
ax.set_xticks(np.arange(90,181,30), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(-50,11,10), crs=ccrs.PlateCarree())
#plt.contour(lonproj, latproj, PP[3,:,:], levels=[0.75], colors='0', \
#            label='0.75 ZFL')
#plt.contour(lonproj, latproj, (PP[1,:,:]/PP[0,:,:]), levels=[1.46], colors='r')
#plt.contour(lonproj, latproj, PP_[2,:,:], levels=[PP_[3,0,0]], colors='0.6', \
#            label='F90 vSZ')
#plt.title('PP ratio - TMM based on ZFL with vSZ ratio (significance based on vSZ in grey)', \
#          fontsize=16, y=1.08)
plt.title('PPR from daily SSTA using vSW')
plt.savefig(figfile,bbox_inches='tight', format='png', dpi=300)
plt.show(block=False)


