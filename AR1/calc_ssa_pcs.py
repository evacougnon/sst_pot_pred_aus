'''
    compute a single spectral analysis on PCs time series
    determined after applying and EOF analysis (using cdo) 
    on the SSTa 

    Following the idea in Monselesan et al 2015 and O'Kane et al 2016

    With the help of:
    https://dept.atmos.ucla.edu/tcd/ssa-tutorial-matlab

    Author: Eva C.
    Created: Oct 2018
    Last Modif:


    TESTING phase!
    What about comparing, writing it and using the mcssa toolbox 
    for SSA analysis i python, also include onte Carlo SSA 
    (https://github.com/VSainteuf/mcssa)



############ check the matlab script, more efficient!!!!!


'''

# library
# import libraries
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt

from mcssa import mcssa # to compute SSA analysis
# NB: tested to compare with tutorial: 
# https://dept.atmos.ucla.edu/tcd/ssa-tutorial-matlab
# and script cloned from https://github.com/VSainteuf/mcssa
# all advertised on UCLA: https://dept.atmos.ucla.edu/tcd/matlab-tutorials

#####################################################
# load data
#####################################################
# EOF on the yearly mean and 1x1 deg resolution grid
'''
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/yearly/PC_yearlyONmonthly_1deg.nc'
M=12 #dimension of the desired embedding (window length)
'''
# EOF on the unfiltered (daily) SSTa and (1/4)x(1/4) deg resolution grid
fname_pc = '/v_Munk_Drive/ecougnon/scripts/cdo_tool/PCs_daily_trend_1deg.nc'
M=360


########################
# choose a mode to do the SSA analysis on:
hmode=1
eof_pcs = xr.open_dataset(fname_pc)['SSTa'][:,hmode-1].squeeze()

#######
# computes SSA 
ssa=mcssa.SSA(eof_pcs)

ssa.run_ssa(M)

# reconstructed components (the sum of all of them 
# should be = to the original time series)
RC_all=np.zeros((M,length(eof_pcs)))
for ii in range(1,M+1):
    RC_all[ii-1,:] = ssa.reconstruct([ii])


##########################
# plotting

# plotting the reconstructed component:
plt.figure()
plt.plot(RC_all[0])
plt.show(block=False)


# checking that the summ of the M reconstructed component is 
# equal to the original time series (PC)
plt.figure()
plt.plot(eof_pcs)
plt.plot(np.sum(RC_all,axis=0))
plt.show(block=False)


