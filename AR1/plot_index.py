'''
   calculate the correlation coefficient (Spearman rank coef)
   of a selected leading mode compared to several indeces:
   MEI, NINO34, NINO3, NINO4, DMI, SAM (marshall 2003), and 
   blocking?!
'''

# libriries
import numpy as np
import xarray as xr
import scipy.stats.mstats as st # highly similar to scipy.stats 
				# but adapted for masked arrays
from scipy import signal
from scipy import stats
import pandas as pd
import csv
from matplotlib import pyplot as plt

# useful numbres
MinYear = 1982
MaxYear = 2016
NumYears = MaxYear-MinYear+1

## Load the selected leading mode
figfile = '/home/ecougnon/Desktop/WorkingFigures/vSZ/EOF/corr_coef_1degAus_trend_PC8.png'
fname = '../../ana/PotPred/EOF/eof_Aus_daily_trend_1deg_12modes_monthly.nc'
tim_vec = xr.open_dataset(fname)['time'] 
#data = np.load(fname + '.npz')
#PC = data['PC'].item()

# which mode to compare with
mode = 8
print(mode)
var_PC = xr.open_dataset(fname)['PCs'][-mode,:] #PC['PC1_pred'].copy()

## CLIMATE MODE
# read MEI index (ENSO)
# file_enso is a file with 67 (axis=0) years (starting in 1950) and 
# 12 moonths + 1 year label columns (13, axis=1) 
file_mei = np.genfromtxt('../../data/index/enso/MEI_index.txt', \
                         skip_header=10, skip_footer = 30, delimiter='\t')
str_mei = np.nonzero((file_mei[:,0]>(MinYear-1)) \
                     & (file_mei[:,0]<(MinYear+1)))
mei_monthly = np.empty(NumYears*12)
mei_monthly_std = np.empty(NumYears*12)
k=0
for yy in np.arange(str_mei[0][0],len(file_mei[:,0])):
    for mm in np.arange(1,12+1):
        mei_monthly[k] = file_mei[yy,mm]
        k = k + 1
mei_monthly_std = (mei_monthly-np.nanmean(mei_monthly)) \
                  /np.nanstd(mei_monthly)
mei_monthly = signal.detrend(mei_monthly)
#NINO34
file_nino34 = np.genfromtxt('../../data/index/enso/nino34.txt', \
                            skip_header=1, skip_footer = 6)
file_nino34_a = np.genfromtxt('../../data/index/enso/nino34_anomaly.txt', \
                              skip_header=1, skip_footer = 7)

str_nino34 = np.nonzero((file_nino34[:,0]>(MinYear-1)) \
                        & (file_nino34[:,0]<(MinYear+1)))
nino34_monthly = np.empty(NumYears*12)
nino34_monthly_a = np.empty(NumYears*12)
nino34_monthly_std = np.empty(NumYears*12)
k=0
for yy in np.arange(str_nino34[0][0],len(file_nino34[:,0])):
    for mm in np.arange(1,12+1):
        nino34_monthly[k] = file_nino34[yy,mm]
        nino34_monthly_a[k] = file_nino34_a[yy,mm]
        k = k + 1
nino34_monthly = signal.detrend(nino34_monthly)
# standardise the index
nino34_monthly_std = (nino34_monthly-np.nanmean(nino34_monthly)) \
                     /np.nanstd(nino34_monthly)
# NINO3
file_nino3 = np.genfromtxt('../../data/index/enso/nino3.txt', \
                           skip_header=1, skip_footer = 6)
file_nino3_a = np.genfromtxt('../../data/index/enso/nino3_anomaly.txt', \
                             skip_header=1, skip_footer = 7)

str_nino3 = np.nonzero((file_nino3[:,0]>(MinYear-1)) \
                       & (file_nino3[:,0]<(MinYear+1)))
nino3_monthly = np.empty(NumYears*12)
nino3_monthly_a = np.empty(NumYears*12)
nino3_monthly_std = np.empty(NumYears*12)
k=0
for yy in np.arange(str_nino3[0][0],len(file_nino3[:,0])):
    for mm in np.arange(1,12+1):
        nino3_monthly[k] = file_nino3[yy,mm]
        nino3_monthly_a[k] = file_nino3_a[yy,mm]
        k = k + 1
nino3_monthly = signal.detrend(nino3_monthly)
# standardise the index
nino3_monthly_std = (nino3_monthly-np.nanmean(nino3_monthly)) \
                    /np.nanstd(nino3_monthly)

# NINO4
file_nino4 = np.genfromtxt('../../data/index/enso/nino4.txt', \
                           skip_header=1, skip_footer = 6)
file_nino4_a = np.genfromtxt('../../data/index/enso/nino4_anomaly.txt', \
                             skip_header=1, skip_footer = 7)

str_nino4 = np.nonzero((file_nino4[:,0]>(MinYear-1)) \
                       & (file_nino4[:,0]<(MinYear+1)))
nino4_monthly = np.empty(NumYears*12)
nino4_monthly_a = np.empty(NumYears*12)
nino4_monthly_std = np.empty(NumYears*12)
k=0
for yy in np.arange(str_nino4[0][0],len(file_nino4[:,0])):
    for mm in np.arange(1,12+1):
        nino4_monthly[k] = file_nino4[yy,mm]
        nino4_monthly_a[k] = file_nino4_a[yy,mm]
        k = k + 1
nino4_monthly = signal.detrend(nino4_monthly)
# standardise the index
nino4_monthly_std = (nino4_monthly-np.nanmean(nino4_monthly)) \
                    /np.nanstd(nino4_monthly)

# IOD (DMI index)
file_dmi = np.genfromtxt('../../data/index/DMI_IOD.txt', \
                         skip_header=1, skip_footer = 4)

str_dmi = np.nonzero((file_dmi[:,0]>(MinYear-1)) \
                     & (file_dmi[:,0]<(MinYear+1)))
dmi_monthly = np.empty(NumYears*12)
dmi_monthly_std = np.empty(NumYears*12)
k=0
for yy in np.arange(str_dmi[0][0],len(file_dmi[:,0])):
    for mm in np.arange(1,12+1):
        dmi_monthly[k] = file_dmi[yy,mm]
        k = k + 1
dmi_monthly = signal.detrend(dmi_monthly)
# standardise the index
dmi_monthly_std = (dmi_monthly-np.nanmean(dmi_monthly)) \
                  /np.nanstd(dmi_monthly)
# SAM
file_sam = np.genfromtxt('../../data/index/sam/marshall2003.txt', \
                         skip_header=2, skip_footer = 3)

str_sam = np.nonzero((file_sam[:,0]>(MinYear-1)) \
                     & (file_sam[:,0]<(MinYear+1)))
sam_monthly = np.empty(NumYears*12)
sam_monthly_std = np.empty(NumYears*12)
k=0
for yy in np.arange(str_sam[0][0],len(file_sam[:,0])):
    for mm in np.arange(1,12+1):
        sam_monthly[k] = file_sam[yy,mm]
        k = k + 1
sm_monthly = signal.detrend(sam_monthly)
# standardise the index
sam_monthly_std = (sam_monthly-np.nanmean(sam_monthly)) \
                   /np.nanstd(sam_monthly)

#subtropical ridge tasmanian high
# data from 1982-2012 in the following file
data_strh = np.load('../../data/strh_index.npz')
strh_monthly = data_strh['strh']
strh_lim = len(strh_monthly)

# blocking index centred at 140E
data_bi_140 = np.load('../../data/blocking_index_140.npz')
bi_140_monthly = data_bi_140['BI_140']
bi_140_monthly = signal.detrend(bi_140_monthly)
# blocking index centred at 160E
data_bi_160 = np.load('../../data/blocking_index_160.npz')
bi_160_monthly = data_bi_160['BI_160']
bi_160_monthly = signal.detrend(bi_160_monthly)

# Ningaloo nino index
data_NingN = np.load('../../data/NingalooNino_index.npz')
ning_monthly = data_NingN['NingN_Ind']
ning_monthly = signal.detrend(ning_monthly)

'''
# EAC transport 
data_eac= np.load('../../data/EAC_trp.npz')
eac_time = data_eac['tim_vec']
eac_monthly = data_eac['trp_eac']
eac_str = int(np.array(np.where(tim_vec == eac_time[0])))
eac_end = int(np.array(np.where(tim_vec == eac_time[-1])))
eac_monthly = signal.detrend(eac_monthly)
'''

# EAC transport -- BRAN from Zeya's analysis
# 1994-2016 (Aug)
df = pd.read_csv('../../data/transport_y.csv', header=None)
eac_monthly = (df.iloc[:,0])*110000*np.cos((np.pi/180)*37)*0.1*10**(-6)
eac_monthly = signal.detrend(eac_monthly)
eac_str = int(12*12)
eac_end = int(12*12+len(eac_monthly)-1)

# calculate correlation coefficient using the spearman rank
corr_mei, p_mei = st.spearmanr(var_PC, mei_monthly)
corr_mei_std, p_mei_std = st.spearmanr(var_PC,mei_monthly_std)

corr_nino34, p_nino34 = st.spearmanr(var_PC,nino34_monthly)
corr_nino34_a, p_nino34_a = st.spearmanr(var_PC,nino34_monthly_a)
corr_nino34_std, p_nino34_std = st.spearmanr(var_PC,nino34_monthly_std)

corr_nino3, p_nino3 = st.spearmanr(var_PC,nino3_monthly)
corr_nino3_a, p_nino3_a = st.spearmanr(var_PC,nino3_monthly_a)
corr_nino3_std, p_nino3_std = st.spearmanr(var_PC,nino3_monthly_std)

corr_nino4, p_nino4 = st.spearmanr(var_PC,nino4_monthly)
corr_nino4_a, p_nino4_a = st.spearmanr(var_PC,nino4_monthly_a)
corr_nino4_std, p_nino4_std = st.spearmanr(var_PC,nino4_monthly_std)

corr_dmi, p_dmi = st.spearmanr(var_PC,dmi_monthly)
corr_dmi_std, p_dmi_std = st.spearmanr(var_PC,dmi_monthly_std)

corr_sam, p_sam = st.spearmanr(var_PC,sam_monthly)
corr_sam_std, p_sam_std = st.spearmanr(var_PC,sam_monthly_std)

corr_strh, p_strh = st.spearmanr(var_PC[0:strh_lim],strh_monthly)
#corr_strh_std, p_strh_std = st.spearmanr(var_PC,strh_monthly_std)

corr_bi_140, p_bi_140 = st.spearmanr(var_PC[0:strh_lim],bi_140_monthly)
corr_bi_160, p_bi_160 = st.spearmanr(var_PC[0:strh_lim],bi_160_monthly)

corr_ning, p_ning = st.spearmanr(var_PC, ning_monthly)

corr_eac, p_eac = st.spearmanr(var_PC[eac_str:eac_end+1], eac_monthly)

# calc correlation coefficient with lead time on the indeces
# do the salculation every month up to 2 years lead 
lead_max = 48+1 #24+1 # lead calculation up to 24 months 
corr_mei_lead = np.empty(len(range(0,lead_max)))
corr_mei_lead.fill(np.nan)
corr_mei_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, mei_monthly) # zero lead/lag

corr_nino34_lead = np.empty(len(range(0,lead_max)))
corr_nino34_lead.fill(np.nan)
corr_nino34_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, nino34_monthly) # zero lead/lag
corr_nino34_a_lead = np.empty(len(range(0,lead_max)))
corr_nino34_a_lead.fill(np.nan)
corr_nino34_a_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, nino34_monthly_a) # zero lead/lag

corr_nino3_lead = np.empty(len(range(0,lead_max)))
corr_nino3_lead.fill(np.nan)
corr_nino3_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, nino3_monthly) # zero lead/lag
corr_nino3_a_lead = np.empty(len(range(0,lead_max)))
corr_nino3_a_lead.fill(np.nan)
corr_nino3_a_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, nino3_monthly_a) # zero lead/lag

corr_nino4_lead = np.empty(len(range(0,lead_max)))
corr_nino4_lead.fill(np.nan)
corr_nino4_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, nino4_monthly) # zero lead/lag
corr_nino4_a_lead = np.empty(len(range(0,lead_max)))
corr_nino4_a_lead.fill(np.nan)
corr_nino4_a_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, nino4_monthly_a) # zero lead/lag

corr_dmi_lead = np.empty(len(range(0,lead_max)))
corr_dmi_lead.fill(np.nan)
corr_dmi_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, dmi_monthly) # zero lead/lag

corr_sam_lead = np.empty(len(range(0,lead_max)))
corr_sam_lead.fill(np.nan)
corr_sam_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, sam_monthly) # zero lead/lag

corr_strh_lead = np.empty(len(range(0,lead_max)))
corr_strh_lead.fill(np.nan)
corr_strh_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC[0:strh_lim], strh_monthly) # zero lead/lag

corr_bi_140_lead = np.empty(len(range(0,lead_max)))
corr_bi_140_lead.fill(np.nan)
corr_bi_140_lead[int(lead_max/2)],tmp = st.spearmanr(var_PC[0:strh_lim],bi_140_monthly) #zero lead/lag
corr_bi_160_lead = np.empty(len(range(0,lead_max)))
corr_bi_160_lead.fill(np.nan)
corr_bi_160_lead[int(lead_max/2)],tmp = st.spearmanr(var_PC[0:strh_lim],bi_160_monthly) #zero lead/lag

corr_ning_lead = np.empty(len(range(0,lead_max)))
corr_ning_lead.fill(np.nan)
corr_ning_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC, ning_monthly) # zero lead/lag

corr_eac_lead = np.empty(len(range(0,lead_max)))
corr_eac_lead.fill(np.nan)
corr_eac_lead[int(lead_max/2)], tmp = st.spearmanr(var_PC[eac_str:eac_end+1], \
                                                   eac_monthly) # zero lead/lag


ll=1
for lag in range(int(lead_max/2)-1,0-1,-1):
#    print(lag)
    corr_mei_lead[lag], tmp = st.spearmanr(var_PC[:-ll],mei_monthly[ll:])
    corr_nino34_lead[lag], tmp = st.spearmanr(var_PC[:-ll], \
                                              nino34_monthly[ll:])
    corr_nino34_a_lead[lag], tmp = st.spearmanr(var_PC[:-ll], \
                                                nino34_monthly_a[ll:])
    corr_nino3_lead[lag], tmp = st.spearmanr(var_PC[:-ll], \
                                             nino3_monthly[ll:])
    corr_nino3_a_lead[lag], tmp = st.spearmanr(var_PC[:-ll], \
                                               nino3_monthly_a[ll:])
    corr_nino4_lead[lag], tmp = st.spearmanr(var_PC[:-ll], \
                                             nino4_monthly[ll:])
    corr_nino4_a_lead[lag], tmp = st.spearmanr(var_PC[:-ll], \
                                               nino4_monthly_a[ll:])
    corr_dmi_lead[lag], tmp = st.spearmanr(var_PC[:-ll],dmi_monthly[ll:])
    corr_sam_lead[lag], tmp = st.spearmanr(var_PC[:-ll],sam_monthly[ll:])
    corr_strh_lead[lag], tmp = st.spearmanr(var_PC[:strh_lim-ll], \
                                            strh_monthly[ll:])
    corr_bi_140_lead[lag], tmp = st.spearmanr(var_PC[:strh_lim-ll], \
                                              bi_140_monthly[ll:])
    corr_bi_160_lead[lag], tmp = st.spearmanr(var_PC[:strh_lim-ll], \
                                              bi_160_monthly[ll:])
    corr_ning_lead[lag], tmp = st.spearmanr(var_PC[:-ll],ning_monthly[ll:])
    corr_eac_lead[lag], tmp = st.spearmanr(var_PC[eac_str:eac_end-ll+1], \
                                           eac_monthly[ll:])

    ll = ll+1

ll=1
for lead in range(int(lead_max/2)+1,lead_max): 
#    print(lead)
    corr_mei_lead[lead], tmp = st.spearmanr(var_PC[ll:],mei_monthly[:-ll]) 
    corr_nino34_lead[lead], tmp = st.spearmanr(var_PC[ll:], \
                                               nino34_monthly[:-ll])
    corr_nino34_a_lead[lead], tmp = st.spearmanr(var_PC[ll:], \
                                                 nino34_monthly_a[:-ll])
    corr_nino3_lead[lead], tmp = st.spearmanr(var_PC[ll:], \
                                              nino3_monthly[:-ll])
    corr_nino3_a_lead[lead], tmp = st.spearmanr(var_PC[ll:], \
                                                nino3_monthly_a[:-ll])
    corr_nino4_lead[lead], tmp = st.spearmanr(var_PC[ll:], \
                                              nino4_monthly[:-ll])
    corr_nino4_a_lead[lead], tmp = st.spearmanr(var_PC[ll:], \
                                                nino4_monthly_a[:-ll])
    corr_dmi_lead[lead], tmp = st.spearmanr(var_PC[ll:],dmi_monthly[:-ll])
    corr_sam_lead[lead], tmp = st.spearmanr(var_PC[ll:],sam_monthly[:-ll])
    corr_strh_lead[lead], tmp = st.spearmanr(var_PC[ll:strh_lim], \
                                             strh_monthly[:-ll])
    corr_bi_140_lead[lead], tmp = st.spearmanr(var_PC[ll:strh_lim], \
                                               bi_140_monthly[:-ll])
    corr_bi_160_lead[lead], tmp = st.spearmanr(var_PC[ll:strh_lim], \
                                               bi_160_monthly[:-ll])
    corr_ning_lead[lead], tmp = st.spearmanr(var_PC[ll:],ning_monthly[:-ll])
    corr_eac_lead[lead], tmp = st.spearmanr(var_PC[eac_str+ll:eac_end+1], \
                                            eac_monthly[:-ll])

    ll = ll+1
 
print('mei: ' + repr(corr_mei) + '' )
#print('mei std: ' + repr(corr_mei_std) + '')

print('nino34: '+ repr(corr_nino34) +'')
print('nino34 a: '+ repr(corr_nino34_a) +'')
#print('nino34 std: '+ repr(corr_nino34_std) +'')

print('nino3: '+ repr(corr_nino3) +'')
print('nino3 a: '+ repr(corr_nino3_a) +'')
#print('nino3 std: '+ repr(corr_nino3_std) +'')

print('nino4: '+ repr(corr_nino4) +'')
print('nino4 a: '+ repr(corr_nino4_a) +'')
#print('nino4 std: '+ repr(corr_nino4_std) +'')

print('dmi: '+ repr(corr_dmi) +'')
#print('dmi std: '+ repr(corr_dmi_std) +'')

print('sam: '+ repr(corr_sam) +'')
#print('sam std: '+ repr(corr_sam_std) +'')

print('strh: '+ repr(corr_strh) +'')

print('BI at 140E: '+ repr(corr_bi_140) +'')
print('BI at 160E: '+ repr(corr_bi_160) +'')

print('Ningaloo Nino: '+ repr(corr_ning) +'')

print('EAC transport: '+ repr(corr_eac) +'')

## significance:
N = len(var_PC)
t_interval_max = np.empty(lead_max)
t_interval_max.fill(np.nan)
t_interval_min = np.empty(lead_max)
t_interval_min.fill(np.nan)
t_interval_min[int(lead_max/2)], t_interval_max[int(lead_max/2)] = \
stats.t.interval(0.95,N-1)
for df in range(1,int(lead_max/2)+1):
    N_ = N-df
    t_interval_min[int(lead_max/2)+df], t_interval_max[int(lead_max/2)+df] = \
    stats.t.interval(0.95,N_-1)
    t_interval_min[int(lead_max/2)-df], t_interval_max[int(lead_max/2)-df] = \
    stats.t.interval(0.95,N_-1)
## testing with rs: corr coeff
# rs=0.0985
# rs*np.sqrt((N-df-2)/(1-rs**2))
# 1.9647, while t_interval for the same df = 1.966 (for actually 
# alomost all the df considered here (1.9647-1.9660)
## ploting the corr coef funcion of  the lead time
x_axis = np.arange(-24,24+1)

ax=plt.figure(figsize=(15,6))
plt.plot(x_axis, corr_mei_lead,'c') #'b')
plt.plot(x_axis, corr_nino34_lead,'g')
#plt.plot(x_axis, corr_nino34_a_lead,'g')
##plt.plot(x_axis, corr_nino3_lead,'y')
#plt.plot(x_axis, corr_nino3_a_lead,'y')
##plt.plot(x_axis, corr_nino4_lead,'c')
#plt.plot(x_axis, corr_nino4_a_lead,'c')
plt.plot(x_axis, corr_dmi_lead,'m')
plt.plot(x_axis, corr_sam_lead,'r')
#plt.plot(x_axis, corr_strh_lead,'k')
##plt.plot(corr_bi_140_lead)
#plt.plot(x_axis, corr_bi_160_lead, color='0.75')
plt.plot(x_axis, corr_ning_lead,'darkorange')
plt.plot(x_axis, corr_eac_lead,'k') #'maroon')
##plt.legend(['mei','nino34','nino34a','nino3','nino3a','nino4','nino4a', \
##            'dmi','sam','strh','BI 160'],loc=4)
plt.legend(['mei','nino34','dmi','sam', \
            'Ningaloo Nino', 'EAC trp'],loc=4)
#plt.legend(['Multivariate ENSO Index (MEI)','Dipole Mode Index (DMI)', \
#            'Southern Annular Mode (SAM)','Ningaloo Nino Index (NNI)', \
#            'EAC transport'],loc=2, fontsize=18)
plt.title('correlation coefficient with lead time for PC' + str(int(mode)) + '', \
          fontsize=16)
plt.xlabel('lead month', fontsize=20)
plt.ylabel('Spearman rank coefficient',fontsize=20)
plt.plot([0, 24], [0.0985, 0.0985],'--k')
plt.plot([0, 24], [-0.0985, -0.0985],'--k')
plt.xlim([0, 24])
plt.ylim([-0.6, 0.6])
plt.grid()
plt.tick_params(labelsize=18)

plt.savefig(figfile, bbox_inches='tight', format='png', dpi=300)

#plt.show(block=False)

##
'''
# plotting
plt.figure()
#plt.plot(enso_monthly)
plt.plot(enso_monthly_a)
plt.plot(enso_monthly_std)
plt.legend(['anomalies', 'standardized'])
plt.title('NINO34 index from NOAA')
plt.grid()
plt.show()
'''


