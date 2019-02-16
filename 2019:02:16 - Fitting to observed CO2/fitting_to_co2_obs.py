import numpy as np
import scipy as sp
from scipy.optimize import minimize
import string
import math
import sys
import pandas as pd

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

from fiveEqSCM_multiparam import *

if __name__=='__main__':
	mode = sys.argv[1]

	if mode == 'tune_GCP':

		do_fit = True

		gcp_ems = np.genfromtxt('./data/GCP_ems_2018.txt', delimiter='\t', skip_header=1)

		gcp_ems_plusHistorical = np.genfromtxt('./data/GCP_ems_1850_to_2017.txt', delimiter='\t')

		noaa_co2_annual_gl = np.genfromtxt('./data/noaa_co2_conc.txt', delimiter='\t', skip_header=57)

		RCP85_emissions = np.genfromtxt('./RCP_data/RCP85_EMISSIONS.csv', skip_header=37, delimiter = ',')

		RCP85_rf = np.genfromtxt('./RCP_data/RCP85_MIDYEAR_RADFORCING.csv', skip_header=59, delimiter = ',')

		piers_updated_rf = np.genfromtxt('./data/Annualforcings_Mar2014_GHGrevised.txt', skip_header=4, delimiter = '\t')

		other_rf = piers_updated_rf[1765-1750:,13] - piers_updated_rf[1765-1750:,1]
		other_rf = np.append(other_rf, other_rf[-1])

		# import the HadCRUT4 temperature observations as input to CO2 fit for main CO2 fit?
		# for CO2 fit to MAGICC, use the etminan forcing value, and the temperature response from CO2 emissions plus other RF?

		CO2_ems = np.zeros(int(gcp_ems_plusHistorical[-1,0] - RCP85_emissions[0,0] + 1))
		CO2_ems[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0])] = (RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2]) * (gcp_ems_plusHistorical[0,1] + gcp_ems_plusHistorical[0,2]) / (RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2])
		CO2_ems[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]):] = gcp_ems_plusHistorical[:,1] + gcp_ems_plusHistorical[:,2]

		# other_rf = RCP85_rf[:2018,1] - RCP85_rf[:2018,8]

		C, RF, T = multiscen_oxfair(emissions=CO2_ems,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[32.40,0.019,4.165,0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = 1.6,
									ecs = 2.75,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)
		
		C_rcp, RF_rcp, T_rcp = multiscen_oxfair(emissions=RCP85_emissions[:int(gcp_ems[-1,0] - RCP85_emissions[0,0])+1,1] + RCP85_emissions[:int(gcp_ems[-1,0] - RCP85_emissions[0,0])+1,2],
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[32.40,0.019,4.165,0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = 1.6,
									ecs = 2.75,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)


		if do_fit == True:
			def noaa_diff_co2(x):
				C_vals = multiscen_oxfair(emissions=CO2_ems,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[x[0],x[1],x[2],0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = 1.6,
									ecs = 2.75,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)[0]

				diff = np.sum((C_vals[0,0,1980-1765:2018-1765] - noaa_co2_annual_gl[:,1])**2)
				return diff
			x_val = minimize(noaa_diff_co2, (32.4,0.019,4.165), bounds=((20, 50), (0, 1), (0,15)), method='TNC', options={})

			print(x_val)


		r0_fit = x_val.x

		C_fit, RF_fit, T_fit = multiscen_oxfair(emissions=CO2_ems,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[r0_fit[0],r0_fit[1],r0_fit[2],0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = 1.6,
									ecs = 2.75,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)



		fig1, ax1 = plt.subplots(2,2,figsize=(16,15))

		ax1[0,0].plot(RCP85_emissions[:2019-1765,0], RCP85_emissions[:2019-1765,1] + RCP85_emissions[:2019-1765,2], color = 'red')
		ax1[0,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), CO2_ems, color='green')
		ax1[0,0].plot(gcp_ems[:,0], gcp_ems[:,1] + gcp_ems[:,2], color='blue')
		ax1[0,0].plot(gcp_ems_plusHistorical[:,0], gcp_ems_plusHistorical[:,1] + gcp_ems_plusHistorical[:,2], color='blue')
		ax1[0,0].set_xlabel('Year')
		ax1[0,0].set_ylabel('Annual CO$_2$ Emissions (GtC/yr)')

		ax1[0,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), C_rcp[0,0,:], color='red')
		ax1[0,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), C[0,0,:], color='green')
		ax1[0,1].plot(noaa_co2_annual_gl[:,0], noaa_co2_annual_gl[:,1], color='blue')
		ax1[0,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), C_fit[0,0,:], color='orange', linestyle='--')
		ax1[0,1].set_xlabel('Year')
		ax1[0,1].set_ylabel('CO$_2$ concentration (ppmv)')

		ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), RF_rcp[0,:], color='red')
		ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), RF[0,:], color='green')
		ax1[1,0].set_xlabel('Year')
		ax1[1,0].set_ylabel('CO$_2$ radiative forcing (Wm$^{-2}$)')

		ax1[1,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), T_rcp[0,:] - np.mean(T_rcp[0,1850-1765:1901-1765]), color='red')
		ax1[1,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), T[0,:] - np.mean(T[0,1850-1765:1901-1765]), color='green')
		ax1[1,1].set_xlabel('Year')
		ax1[1,1].set_ylabel('Temperature relative to 1850-1900 (K)')

		ax1[1,1].axhline(y = 1.0, linestyle=':', color='black', alpha=0.4)
		ax1[1,1].axhline(y = 0, linestyle='-', color='black', alpha=0.4)

		plt.show()

