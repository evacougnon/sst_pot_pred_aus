'''

plot maps of the persistence of a forecast using SSTs around 
Australia, after calculating tau from the OTEs_ts_AR1.py
that takes about 45 minutes to run for the Aus region (-55:10N 
and 90:180E)

'''
# import library
import numpy as np

#import matplotlib
#matplotlib.use("TkAgg") # used (once) to allow figures display
from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm


#
# Load data
#
fname = '/home/ecougnon/ana/AR1/SST_annual_extremes_Aus'
data = np.load(fname+'.npz')
lon_map = data['lon_map']
lat_map = data['lat_map']
years = data['years']
SST = data['SST'].item()

fname_tau = '/home/ecougnon/ana/AR1/OTEs_tau_Aus'
data_tau = np.load(fname_tau+'.npz')
a = data_tau['a'].item() # auto-regressive parameter
tau = data_tau['tau'].item() # persitence of the forecast -- time scale
sig = data_tau['sig'].item() 

figfile = '/home/ecougnon/ana/AR1/'

# plotting parameters
# Maps
domain = [-55, 90, 10, 180] # [-80, 20, 85, 380] #[-55, 90, 10, 180] 
domain_draw = [-50, 90, 10, 180] #[-80, 20, 85, 380] #[-55, 90, 10, 180] 
dlat = 10 #30 #10
dlon = 30 #90 #30
llon, llat = np.meshgrid(lon_map, lat_map)
bg_col = '0.6'
cont_col = '1.0'


plt.figure(figsize=(15,17))
for key in SST.keys():
    plt.clf()
    plt.subplot(3,1,1, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#    proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                       labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                       labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    plt.contourf(lonproj, latproj, a[key], levels=np.arange(-1,1+0.2,0.2), \
                 cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    plt.contour(lonproj, latproj, a[key], levels=[0], colors='k')
    plt.title(key + ': AR1 coefficient')
    plt.subplot(3,1,2, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#    proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                       labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                       labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    plt.contourf(lonproj, latproj, tau[key], levels=np.append(np.arange(0,2+0.25,0.25),3), \
                 cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    H = plt.colorbar()
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], \
                 levels=[0., 0.5, 1.0], colors='none')
    H.set_label(r'[years]')
    plt.title(key + r': AR1 timescale')
    plt.subplot(3,1,3, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#    proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                       labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                       labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    plt.contourf(lonproj, latproj, sig[key], levels=np.append(np.arange(0,1.35+0.15,0.15),1.5), \
                 cmap=plt.cm.afmhot_r)
    plt.clim(0,1.3)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.title(key + ': AR1 error std. dev.')
    
#    plt.savefig(figfile + 'OTEs_AR1_' + key + '.png', bbox_inches='tight', \
#                dpi=300) #pad_inches=0.5)

# Annual plots
fig = plt.figure(figsize=(12,5)) #(22,4))
titles = ['Annual Mean', 'Annual P90', 'Annual P10']
cnt = 0
for key in ['TMM', 'TA90', 'TA10']:
    cnt += 1
    AX = plt.subplot(1,3,cnt, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#    proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                       labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                       labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    H = plt.contourf(lonproj, latproj, tau[key], \
                     levels=np.append(np.arange(0,2+0.25,0.25),3), \
                     cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], \
                 levels=[0., 0.5, 1.0], colors='none')
    plt.title(titles[cnt-1])

AXPOS = AX.get_position()
#CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
CAX = fig.add_axes([AXPOS.x1+0.015, 0.27, 0.01, 0.45]) 
HB = plt.colorbar(H,CAX, orientation='vertical') 
HB.set_label(r'[years]')

plt.savefig(figfile + 'OTEs_AR1_Annual.eps', bbox_inches='tight', \
            dpi=300, format='eps') #pad_inches=0.5, dpi=150)

# Seasonal 90th percentile
fig = plt.figure(figsize=(11,7))
titles = ['Summer (P90)', 'Autumn (P90)', 'Winter (P90)', 'Spring (P90)']
cnt = 0
for key in ['TSu90', 'TAu90', 'TWi90', 'TSp90']:
    cnt += 1
    AX = plt.subplot(2,2,cnt, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#    proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                       labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                       labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    H = plt.contourf(lonproj, latproj, tau[key], 
                     levels=np.append(np.arange(0,2+0.25,0.25),3), \
                     cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], \
                 levels=[0., 0.5, 1.0], colors='none')
    plt.title(titles[cnt-1])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[years]')

#plt.savefig(figfile + 'OTEs_AR1_Seasonal_P90.png', bbox_inches='tight', \
#            dpi=300) #pad_inches=0.5, dpi=150)

# Seasonal 10th percentile
fig = plt.figure(figsize=(11,7))
titles = ['Summer (P10)', 'Autumn (P10)', 'Winter (P10)', 'Spring (P10)']
cnt = 0
for key in ['TSu10', 'TAu10', 'TWi10', 'TSp10']:
    cnt += 1
    AX = plt.subplot(2,2,cnt, axisbg=bg_col)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#    proj = bm.Basemap(projection='cyl', llcrnrlat=domain[0], llcrnrlon=domain[1], \
#                      urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                       labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                       labels=[False,False,False,True])
    lonproj, latproj = proj(llon, llat)
    H = plt.contourf(lonproj, latproj, tau[key], \
                     levels=np.append(np.arange(0,2+0.25,0.25),3), \
                     cmap=plt.cm.inferno_r)
    plt.clim(0,2)
    plt.contourf(lonproj, latproj, (a[key]<0).astype(float), hatches=['', '.'], \
                 levels=[0., 0.5, 1.0], colors='none')
    plt.title(titles[cnt-1])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[years]')

#plt.savefig(figfile + 'OTEs_AR1_Seasonal_P10.png', bbox_inches='tight', \
#            dpi=300) #pad_inches=0.5, dpi=150)


