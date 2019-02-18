## Code to generate the scenarios -- stitch historical shapes on to future shapes for CO2 emissions, non-CO2 forcing

# Written by Stuart Jenkins (stuart.jenkins@wadham.ox.ac.uk)

# ---------------------------------------------------------------------------------------------------

# Imports 
# -------------------------------------------------
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.special import erf, erfinv
from scipy.optimize import root
from scipy.stats import norm, beta
from statsmodels.api import OLS
import statsmodels.tools.tools
import datetime
import os
import h5py
import copy
from pandas import DataFrame
# -------------------------------------------

# Import the FaIR climate-carbon-cycle model from python file
from fair_scm import *
from forcing_imports import express_as_anom
# -------------------------------------------------

# Definitions
# -------------------------------------------------
a=3.74/np.log(2.0)
c0=278.
# set baseline period for FaIR model runs
base_low=1861.
base_high=1880.
# Set FaIR model temperature boxes parameters
d=np.array([409.5,8.4])
years = np.arange(1765.0,2501.0)
# -------------------------------------------------






# -----------------------------------------------
# Produces scenario shape depending on inputs
#      co2_times - array of times for CO2 fractions
#      co2_fracs - array of fractions of 2020 annual CO2 emissions (i.e. 0.0 = no emissions at that year, 1.0 = 2020 annual emissions rate at that year)
#      tcrecs_index - which tcrecs values should we sample (1 = mean, 0/2 = likely range, 3/4 = central tercile)
#      rf_targs - array of non-CO2 rf targets, if above 0.7 sets constant non-CO2 rf after 2030, if below 0.5 sets declining non-CO2 rf after 2030.
# -----------------------------------------------
def linear_decline_tozero(TCR, ECS, sf_gauss, sf_aero, rf_scens, rf_cont, rf_nat, iirf_100_preind, s_temp, s_cu, slcp_mode='parab',rf_targs=[0.65],
                               co2_times=np.array([2020,2035,2045,2060,2075]),
                               co2_fracs=1.0-np.array([0,0.025*15,0.025*25,0.025*40,1.0]), tcrecs_index=[1],
                               temp_zero={}, llems_zero={}, slcp_zero={}):

    slcp_targs=rf_targs

    #Load scenario pathways
    # Gives shape of non-CO2 forcing after 2030
    f_dat = np.genfromtxt('./Data/F_oth.txt',skip_header=1)
    f_dat_years = f_dat[:,0]

    for name in ['RCP8.5']:

        #Create the arrays for the idealised scenarios
        if temp_zero == {}:
            temp_zero[name] = np.zeros((len(slcp_targs),len(TCR),len(years)))
            llems_zero[name] = np.zeros((len(slcp_targs),len(TCR),len(years)))
            slcp_zero[name] = np.zeros((len(slcp_targs),len(TCR),len(years)))

        for z in range(0,len(slcp_targs)):

            #Examine the central forcing scenario for scenario creation
            cent_forc = sf_gauss[1] *rf_scens[name]['rf_gauss']+\
                        rf_scens[name]['rf_bc'] + sf_aero[1]*rf_scens[name]['rf_aero'] +\
                        rf_cont

            #Set start date for divergence from RCP8.5 trends     
            # Pathways follow RCP8.5 trajectories until 2020
            parab_start = 2020.

            #Substitute the scenarios between 2020 and 2030 for a high forcing target
            if slcp_targs[z] > 0.7:
                cent_forc[np.logical_and(years>=2020,years<=2030)] = cent_forc[years==2020] + f_dat[np.logical_and(f_dat_years>=2020,f_dat_years<=2030),2] - f_dat[f_dat_years==2020,2][0]
                #Set flat forcing after 2030
                cent_forc[years>=2030] = cent_forc[years==2030] 

            #Substitute the scenarios pathway for the low forcing target
            if slcp_targs[z] < 0.5:

                cent_forc[np.logical_and(years>=2020,years<=2100)] = cent_forc[years==2020] + f_dat[f_dat_years>=2020,2] - f_dat[f_dat_years==2020,2][0]
                cent_forc[years>2100] = cent_forc[years==2100]

            for i in tcrecs_index:
                
                # set non-CO2 forcing (including scaling aerosol and gaussian RF components)
                forc = sf_gauss[i] *rf_scens[name]['rf_gauss']+\
                        rf_scens[name]['rf_bc'] + sf_aero[i]*rf_scens[name]['rf_aero'] +\
                        rf_cont

                #Back out the emissions for CO2  
                tb, comb_ems = fair_scm_emsback(rf_scens[name]['rf_co2'],other_rf= forc+rf_nat,TCR=float(TCR[i]),
                                                ECS=float(ECS[i]),a=a,rT=float(s_temp[i]),
                                                r0=float(iirf_100_preind[i]),rC = float(s_cu[i]),
                                                C_0=c0,d1=d[1],d2=d[0])
                
                # set the emissions between 2017 and 2020 to rise in a straight line
                for pos_i in range(1,4):
                    comb_ems[2017+pos_i-1765] = comb_ems[2017+pos_i-1766] + 0.04

                #Get the 2020 CO2 emissions value 
                ems_2020 = comb_ems[years==2020]

                #Create the future CO2 emissions scenario
                #Linearly interpolate between the co2_frac and co2_times to set future emissions trajectory
                llems_zero[name][z][i] = (comb_ems).copy()
                llems_zero[name][z][i][years>=2020.] = ems_2020 *np.interp(years[years>=2020],co2_times,co2_fracs)
                llems_zero[name][z][i][years>=co2_times[-1]] = llems_zero[name][z][i][years==co2_times[-1]]

                # set scenario after 2020 depending on the chosen rf_targs value (below 0.5 sets lower forcing target, above 0.7 sets constant after 2030)
                forc[years>=2020] = forc[years==2020.] + (cent_forc[years>=2020]-cent_forc[years==2020])
                
                #Integrate the scenarios to get temperature response
                concs_i, temps_i  = fair_scm(emissions=llems_zero[name][z][i],other_rf=forc + rf_nat,
                                            tcrecs=np.array([float(TCR[i]),float(ECS[i])]),F_2x=a*np.log(2),
                                            rT=float(s_temp[i]),r0=float(iirf_100_preind[i]),
                                            rC = float(s_cu[i]),C_0=c0,d=d)

                temp_zero[name][z][i] = express_as_anom(temps_i,years)
                slcp_zero[name][z][i] = forc

    return temp_zero, llems_zero, slcp_zero

# -----------------------------------------------

