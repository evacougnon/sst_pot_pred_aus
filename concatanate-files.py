'''
	concatanate npz files

'''
import numpy as np
import xarray as xr
import pandas as pd


outfile = '/home/ecougnon/ana/OTEs_NumDays_month_Aus_19822016.nc'
fname1 = '/home/ecougnon/ana/OTEs_p10_NumDays_month_Aus_test'
#fname2 = '/home/ecougnon/ana/OTEs_p10_NumDays_month_Aus_20002016.nc'
fname3 = '/home/ecougnon/ana/OTEs_p90_NumDays_month_Aus_test'
#fname4 = '/home/ecougnon/ana/OTEs_p90_NumDays_month_Aus_20002016.nc'

month_vec = pd.date_range('1982-01-01','2016-12-31',name='time',freq='M')

data1 = np.load(fname1 + '.npz')
#data2 = np.load(fname2 + '.npz')
data3 = np.load(fname3 + '.npz')
#data4 = np.load(fname4 + '.npz')
lat = data1['lat']
lon = data1['lon']
#month_vec1 = data1['month_vec']
#month_vec2 = data2['month_vec']
Pdays_p10_part1 = data1['Pdays_p10']
#Pdays_p10_part2 = data2['Pdays_p10']
Pdays_p90_part1 = data3['Pdays_p90']
#Pdays_p90_part2 = data4['Pdays_p90']

X = len(lon)
Y = len(lat)
tim = len(month_vec) 
NumDays_month = xr.Dataset({'Pdays_p10':(('time','lat','lon'),np.zeros((tim, Y, X))), \
                            'Pdays_p90':(('time','lat','lon'),np.zeros((tim, Y, X)))}, \
                           {'time': month_vec, 'lat': lat, 'lon': lon})

NumDays_month['Pdays_p10'][:,:,:] = Pdays_p10_part1.copy()
#NumDays_month['Pdays_p10'][len(month_vec1):,:,:] = Pdays_p10_part2.copy()
NumDays_month['Pdays_p90'][:,:,:] = Pdays_p90_part1.copy()
#NumDays_month['Pdays_p90'][len(month_vec1):,:,:] = Pdays_p90_part2.copy()

NumDays_month.to_netcdf(outfile)



