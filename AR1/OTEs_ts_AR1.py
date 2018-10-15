'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to 
  select SST time series around the globe

'''

import numpy as np
from scipy import signal
from scipy import io

import ecoliver as ecj

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

from trendSimAR1 import ar1fit

#
# Load data and make plots
#

outfile = '/home/ecoliver/Desktop/data/NESP/SST_annual_extremes_Aus'

data = np.load(outfile+'.npz')
lon_map = data['lon_map']
lat_map = data['lat_map']
years = data['years']
SST = data['SST'].item()

# Re-map to run 20E to 380E
i_20E = np.where(lon_map>20)[0][0]
lon_map = np.append(lon_map[i_20E:], lon_map[:i_20E]+360)
for key in SST.keys():
    SST[key] = np.append(SST[key][:,i_20E:,:], SST[key][:,:i_20E,:], axis=1)

# Calculate AR1 fit to data
a = {}
tau = {}
tauD = {}
sig = {}
for key in SST.keys():
    a[key] = np.nan*np.zeros(SST[key].shape[0:2])
    tau[key] = np.nan*np.zeros(SST[key].shape[0:2])
    sig[key] = np.nan*np.zeros(SST[key].shape[0:2])
    tauD[key] = np.nan*np.zeros(SST[key].shape[0:2])

for i in range(len(lon_map)):
    print i+1, len(lon_map)
    for j in range(len(lat_map)):
        for key in SST.keys():
            if np.isnan(np.sum(SST[key][j,i,:])) + (np.sum(SST[key][j,i,:])==0):
                continue
            else:
                a[key][j,i], c, sig[key][j,i] = ar1fit(signal.detrend(SST[key][j,i,:]))
                tau[key][j,i] = -1/np.log(a[key][j,i])
                tauD[key][j,i] = (1 + np.abs(a[key][j,i]))/(1 - np.abs(a[key][j,i]))

# Replace tau for a<0 with 0
for key in SST.keys():
    tau[key][a[key]<0] = np.nan #0

# Maps

domain = [-65, 20, 70, 380]
domain_draw = [-60, 60, 60, 380]
dlat = 30
dlon = 60

domain = [-50, 90, 0, 180]
domain_draw = [-50, 90, 0, 180]
dlat = 20
dlon = 30
llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'

plt.figure(figsize=(15,17))
for key in SST.keys():
    plt.clf()
    plt.subplot(3,1,1, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    plt.contourf(lonproj, latproj, a[key], levels=np.arange(-1,1+0.2,0.2), cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    plt.contour(lonproj, latproj, a[key], levels=[0], colors='k')
    plt.title(key + ': AR1 coefficient')
    plt.subplot(3,1,2, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    plt.contourf(lonproj, latproj, tau[key], levels=np.append(np.arange(0, 2+0.25, 0.25), 3), cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    H = plt.colorbar()
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], levels=[0., 0.5, 1.0], colors='none')
    H.set_label(r'[years]')
    plt.title(key + r': AR1 timescale')
    plt.subplot(3,1,3, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    plt.contourf(lonproj, latproj, sig[key], levels=np.append(np.arange(0, 1.35+0.15, 0.15), 1.5), cmap=plt.cm.afmhot_r)
    plt.clim(0,1.3)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.title(key + ': AR1 error std. dev.')
    #
    plt.savefig('/home/ecoliver/Desktop/data/NESP/AR1/OTEs_AR1_' + key + '.png', bbox_inches='tight', pad_inches=0.5)

# Annual plots
fig = plt.figure(figsize=(22,4))
titles = ['Annual Mean', 'Annual P90', 'Annual P10']
cnt = 0
for key in ['TMM', 'TA90', 'TA10']:
    cnt += 1
    AX = plt.subplot(1,3,cnt, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    H = plt.contourf(lonproj, latproj, tau[key], levels=np.append(np.arange(0, 2+0.25, 0.25), 3), cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], levels=[0., 0.5, 1.0], colors='none')
    plt.title(titles[cnt-1])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[years]')

# plt.savefig('/home/ecoliver/Desktop/data/NESP/AR1/OTEs_AR1_Annual.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Seasonal 90th percentile
fig = plt.figure(figsize=(11,7))
titles = ['Summer (P90)', 'Autumn (P90)', 'Winter (P90)', 'Spring (P90)']
cnt = 0
for key in ['TSu90', 'TAu90', 'TWi90', 'TSp90']:
    cnt += 1
    AX = plt.subplot(2,2,cnt, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    H = plt.contourf(lonproj, latproj, tau[key], levels=np.append(np.arange(0, 2+0.25, 0.25), 3), cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], levels=[0., 0.5, 1.0], colors='none')
    plt.title(titles[cnt-1])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[years]')

# plt.savefig('/home/ecoliver/Desktop/data/NESP/AR1/OTEs_AR1_Seasonal_P90.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Seasonal 10th percentile
fig = plt.figure(figsize=(11,7))
titles = ['Summer (P10)', 'Autumn (P10)', 'Winter (P10)', 'Spring (P10)']
cnt = 0
for key in ['TSu10', 'TAu10', 'TWi10', 'TSp10']:
    cnt += 1
    AX = plt.subplot(2,2,cnt, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    H = plt.contourf(lonproj, latproj, tau[key], levels=np.append(np.arange(0, 2+0.25, 0.25), 3), cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], levels=[0., 0.5, 1.0], colors='none')
    plt.title(titles[cnt-1])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[years]')

# plt.savefig('/home/ecoliver/Desktop/data/NESP/AR1/OTEs_AR1_Seasonal_P10.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Some example time series
ii = {}
jj = {}
ii['Ningaloo'] = ecj.find_nearest(lon_map, 112.65)[1]
jj['Ningaloo'] = ecj.find_nearest(lat_map, -21.83)[1]
ii['East. Tas.'] = ecj.find_nearest(lon_map, 148.58)[1]
jj['East. Tas.'] = ecj.find_nearest(lat_map, -43.03)[1]
ii['North. GBR'] = ecj.find_nearest(lon_map, 143.63)[1]
jj['North. GBR'] = ecj.find_nearest(lat_map, -11.21)[1]

plt.figure(figsize=(19,11))
plt.clf()
cnt = 0
for key in ii.keys():
    cnt += 1
    plt.subplot(3,3,cnt)
    plt.plot(years, SST['TMM'][jj[key],ii[key],:], 'k-o')
    plt.plot(years, SST['TA90'][jj[key],ii[key],:], 'r-o')
    plt.plot(years, SST['TA10'][jj[key],ii[key],:], 'b-o')
    plt.grid()
    plt.title(key + ': Annual')
    if key == 'East. Tas.':
        plt.legend(['Mean', 'P90', 'P10'], fontsize=10, loc='upper left')
        plt.ylim(plt.ylim()[0], plt.ylim()[1]+1)
    plt.subplot(3,3,cnt+3)
    plt.plot(years, SST['TWi90'][jj[key],ii[key],:], 'b-o')
    plt.plot(years, SST['TSp90'][jj[key],ii[key],:], 'g-o')
    plt.plot(years, SST['TSu90'][jj[key],ii[key],:], 'r-o')
    plt.plot(years, SST['TAu90'][jj[key],ii[key],:], '-o', color=(1,0.5,0))
    plt.grid()
    plt.title('Seasonal P90')
    if key == 'East. Tas.':
        plt.legend(['Win', 'Spr', 'Sum', 'Aut'], fontsize=10, loc='upper left')
        plt.ylim(plt.ylim()[0], plt.ylim()[1]+1)
    plt.subplot(3,3,cnt+6)
    plt.plot(years, SST['TWi10'][jj[key],ii[key],:], 'b-o')
    plt.plot(years, SST['TSp10'][jj[key],ii[key],:], 'g-o')
    plt.plot(years, SST['TSu10'][jj[key],ii[key],:], 'r-o')
    plt.plot(years, SST['TAu10'][jj[key],ii[key],:], '-o', color=(1,0.5,0))
    plt.grid()
    plt.title('Seasonal P10')
    if key == 'East. Tas.':
        plt.legend(['Win', 'Spr', 'Sum', 'Aut'], fontsize=10, loc='upper left')
        plt.ylim(plt.ylim()[0], plt.ylim()[1]+1)

# plt.savefig('/home/ecoliver/Desktop/data/NESP/AR1/OTEs_AR1_ts.png', bbox_inches='tight', pad_inches=0.5)


