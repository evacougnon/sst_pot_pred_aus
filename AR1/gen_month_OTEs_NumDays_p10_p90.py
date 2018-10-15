'''
	generate a file with the number of days for each month above
	or below the TP90 or TP10 calculated from the entire time series
	of daily SSTa
'''

# import libraries
import numpy as np
import xarray as xr
import pandas as pd
import time as time
import datetime
from dateutil.relativedelta import relativedelta
from calendar import monthrange

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

t_key = time.time()
outfile = '/home/ecougnon/ana/OTEs_p90_NumDays_month_Aus_test' #_19821999'
fname = '/home/ecougnon/ana/SSTa_daily_Aus_shift_time_dim.nc'
fname_p = '/home/ecougnon/ana/OTEs_p10_p90_Aus'
lat = xr.open_dataset(fname)['lat']
lon = xr.open_dataset(fname)['lon']
tim = xr.open_dataset(fname)['time']
tim = tim.sel(time = slice('1982-01-02','2016-12-31')) #'1999-12-31'))
SSTa = xr.open_dataset(fname)['SSTa']
SSTa = SSTa.sel(time=tim, lat=lat, lon=lon)
data = np.load(fname_p + '.npz')
SSTa_p10 = data['SSTa_p90']
#SSTa = SSTa.transpose('time','lat','lon')
#SSTa.to_netcdf('/home/ecougnon/ana/SSTa_daily_Aus_shift_time_dim.nc')
Y = len(lon)
X = len(lat)
# time period
time_vec = pd.date_range('1982-01-01','2016-12-31',name='time')
month_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')
# time index to calc 
MinYear=1982 # 2000
MaxYear=2016 #1999 #2016
NumYears = MaxYear - MinYear+1
MaxNumLeapYear = NumYears//4 + 1 # use only the integer (+1 is used in 
                                 # case the first year is a leap year
NumDays = 365*NumYears + MaxNumLeapYear
dtime = [datetime.date(MinYear,1,1) + datetime.timedelta(days=i) \
         for i in range(NumDays)]
str_id=np.empty(len(month_vec))
end_id=np.empty(len(month_vec))
k=0
for yy in range(MinYear,MaxYear+1):
    for mm in range(1,12+1):
        tmp1 = np.where(dtime == np.datetime64(datetime.date(yy,mm,1)))
        str_id[k] = tmp1[0][0]
        tmp2 = np.where(dtime == np.datetime64(datetime.\
                                               date(yy,mm, \
                                                    monthrange(yy,mm)[1])))
        end_id[k] = tmp2[0][0]
        k = k +1

#SSTa_p10 = np.percentile(SSTa,10,axis=0)

#elapsed_1 = time.time() - t_key
#print('elapsed time to calc p10', elapsed_1)

Pdays_p90 = np.empty((len(month_vec), len(SSTa[0,:,0]), len(SSTa[0,0,:])))
Pdays_p90.fill(np.nan)
k=0
for mm in range(0,len(month_vec)):
#    print('month count: ' + mm)
    for ln in range(0,Y):
        for lt in range(0,X):
#            Pdays_p10[mm,lt,ln] = np.count_nonzero(SSTa[str_id[k]:end_id[k],lt,ln] \
#                                                   < SSTa_p10[lt,ln])
            Pdays_p90[mm,lt,ln] = np.count_nonzero(SSTa[str_id[k]:end_id[k],lt,ln] \
                                                   > SSTa_p10[lt,ln])
    k +=1

elapsed_2 = time.time() - t_key
print('elapsed time to calc day per months', elapsed_2)

np.savez(outfile, Pdays_p90=Pdays_p90, lat=lat, lon=lon, month_vec=month_vec)

'''
# map
domain = [-55, 90, 10, 180] #[-55, 90, 10, 180]
domain_draw = [-50, 90, 10, 180] #[-55, 90, 10, 180]
dlat = 10
dlon = 30
llon, llat = np.meshgrid(lon, lat)
bg_col = '0.6'
cont_col = '1.0'
bin_col = 0.1
bin_bar = 0.5
ax=plt.figure(figsize=(11,11))
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], \
                  urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.fillcontinents(color=(0,0,0), lake_color=(0,0,0), ax=None, \
                    zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), \
                   labels=[True,False,False,False], fontsize=14)
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), \
                   labels=[False,False,False,True], fontsize=14)
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, SSTa_p10, cmap=plt.cm.Reds)

plt.show()
'''



