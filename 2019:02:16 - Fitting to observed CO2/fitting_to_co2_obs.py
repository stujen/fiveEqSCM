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
from Support_codes.temperature_import_and_baseline import *

if __name__=='__main__':
	mode = sys.argv[1]

	if mode == 'tune_GCP':

		do_fit = True

		tcr_val = 1.6
		ecs_val = 2.75

		gcp_ems = np.genfromtxt('./data/GCP_ems_2018.txt', delimiter='\t', skip_header=1)

		gcp_ems_plusHistorical = np.genfromtxt('./data/GCP_ems_1850_to_2017.txt', delimiter='\t')

		noaa_co2_annual_gl = np.genfromtxt('./data/noaa_co2_conc.txt', delimiter='\t', skip_header=57)

		RCP85_emissions = np.genfromtxt('./RCP_data/RCP85_EMISSIONS.csv', skip_header=37, delimiter = ',')

		RCP85_rf = np.genfromtxt('./RCP_data/RCP85_MIDYEAR_RADFORCING.csv', skip_header=59, delimiter = ',')

		piers_updated_rf = np.genfromtxt('./data/Annualforcings_Mar2014_GHGrevised.txt', skip_header=4, delimiter = '\t')

		gmst = temp_import()
		gmst = calc_mean_min_max(gmst)
		annual_temp = np.zeros(168)
		for i in range(0,168):
			annual_temp[i] = np.mean(gmst['Temp-mean'][i*12:i*12 + 12])

		other_rf = 1.0*(piers_updated_rf[1765-1750:,13] - piers_updated_rf[1765-1750:,14]) + 1.0*(piers_updated_rf[1765-1750:,14] - piers_updated_rf[1765-1750:,1])
		other_rf = np.append(other_rf, other_rf[-1])

		# other_rf = (piers_updated_rf[1765-1750:,13] - piers_updated_rf[1765-1750:,1])
		# other_rf = np.append(other_rf, other_rf[-1])

		# import the HadCRUT4 temperature observations as input to CO2 fit for main CO2 fit?
		# for CO2 fit to MAGICC, use the etminan forcing value, and the temperature response from CO2 emissions plus other RF?

		CO2_ems = np.zeros(int(gcp_ems_plusHistorical[-1,0] - RCP85_emissions[0,0] + 1))
		CO2_ems[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0])] = (RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2]) * (gcp_ems_plusHistorical[0,1] + gcp_ems_plusHistorical[0,2]) / (RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2])
		CO2_ems[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]):] = gcp_ems_plusHistorical[:,1] + gcp_ems_plusHistorical[:,2]

		CO2_ems_upper = np.zeros(int(gcp_ems_plusHistorical[-1,0] - RCP85_emissions[0,0] + 1))
		CO2_ems_upper[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]):] = gcp_ems_plusHistorical[:,1]*1.05 + gcp_ems_plusHistorical[:,2] + 0.7
		CO2_ems_upper[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0])] = (RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2]) * (CO2_ems_upper[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0])]) / (RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2])

		CO2_ems_lower = np.zeros(int(gcp_ems_plusHistorical[-1,0] - RCP85_emissions[0,0] + 1))
		CO2_ems_lower[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]):] = gcp_ems_plusHistorical[:,1]*0.95 + gcp_ems_plusHistorical[:,2] - 0.7
		CO2_ems_lower[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0])] = (RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[:int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2]) * (CO2_ems_lower[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0])]) / (RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),1] + RCP85_emissions[int(gcp_ems_plusHistorical[0,0]-RCP85_emissions[0,0]),2])

		other_rf_magicc = RCP85_rf[:2018-1765,1] - RCP85_rf[:2018-1765,8]

		C, RF, T = multiscen_oxfair(emissions=CO2_ems,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[32.40,0.019,4.165,0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = tcr_val,
									ecs = ecs_val,
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
									tcr = tcr_val,
									ecs = ecs_val,
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
									tcr = tcr_val,
									ecs = ecs_val,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)[0]

				diff = np.sum((C_vals[0,0,1980-1765:2018-1765] - noaa_co2_annual_gl[:,1])**2)
				return diff

			def noaa_diff_upper_co2(x):
				C_vals = multiscen_oxfair(emissions=CO2_ems_upper,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[x[0],x[1],x[2],0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = tcr_val,
									ecs = ecs_val,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)[0]

				diff = np.sum((C_vals[0,0,1980-1765:2018-1765] - noaa_co2_annual_gl[:,1])**2)
				return diff

			def noaa_diff_lower_co2(x):
				C_vals = multiscen_oxfair(emissions=CO2_ems_lower,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[x[0],x[1],x[2],0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = tcr_val,
									ecs = ecs_val,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)[0]

				diff = np.sum((C_vals[0,0,1980-1765:2018-1765] - noaa_co2_annual_gl[:,1])**2)
				return diff
			x_val = minimize(noaa_diff_co2, (32.4,0.019,4.165), bounds=((28, 36.6), (0, 0.05), (0,10)), method='TNC', options={'ftol':1e-2, 'xtol':1e-4})

			print(x_val)

			# x_val_upper = minimize(noaa_diff_upper_co2, (32.4,0.019,4.165), bounds=((28, 36.6), (0, 0.05), (0,10)), method='TNC', options={'ftol':1e-2, 'xtol':1e-4})
			# x_val_lower = minimize(noaa_diff_lower_co2, (32.4,0.019,4.165), bounds=((28, 36.6), (0, 0.05), (0,10)), method='TNC', options={'ftol':1e-2, 'xtol':1e-4})

			# print(x_val_lower)
			# print(x_val_upper)

			# print(x_val.x)
			# print(x_val_upper.x)
			# print(x_val_lower.x)


		# r0_fit_upper = np.array([20., 0., 8.77534602])
		# r0_fit_lower = np.array([5.00000000e+01, 9.30489699e-03, 0.00000000e+00])

		r0_fit_upper = np.array([28., 0.0164, 3.602])
		r0_fit_lower = np.array([36.6, 0.0215, 4.706])

		C_fit_upper, RF_fit_upper, T_fit_upper = multiscen_oxfair(emissions=CO2_ems_upper,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[r0_fit_upper[0],r0_fit_upper[1],r0_fit_upper[2],0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = tcr_val,
									ecs = ecs_val,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)

		C_fit_lower, RF_fit_lower, T_fit_lower = multiscen_oxfair(emissions=CO2_ems_lower,
									other_rf = other_rf,
									emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),
									a = np.array([[0.2173,0.2240,0.2824,0.2763]]),
									tau = np.array([[1000000,394.4,36.54,4.304]]),
									r = np.array([[r0_fit_lower[0],r0_fit_lower[1],r0_fit_lower[2],0.0]]),
									PI_C = np.array([278.0]),
									iirf100_max = 97.0,
									f = np.array([[5.78,0.,0.]]),
									tcr = tcr_val,
									ecs = ecs_val,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)




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
									tcr = tcr_val,
									ecs = ecs_val,
									d = np.array([239.0,4.1]),
									q = np.array([0.33,0.41]),
									F_2x = 3.74,
									multigas=False, 
									multiscen=False)



		fig1, ax1 = plt.subplots(2,2,figsize=(10,8))

		# ax1[0,0].plot(RCP85_emissions[:2019-1765,0], RCP85_emissions[:2019-1765,1] + RCP85_emissions[:2019-1765,2], color = 'red')
		ax1[0,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), CO2_ems, color='orange')
		ax1[0,0].plot(gcp_ems_plusHistorical[:,0], gcp_ems_plusHistorical[:,1] + gcp_ems_plusHistorical[:,2], color='black', linestyle='--', alpha=0.5)
		ax1[0,0].set_xlabel('Year')
		ax1[0,0].set_ylabel('Annual CO$_2$ Emissions (GtC/yr)')
		ax1[0,0].text(1870,11,'GCP CO$_2$ ems (1850-2018) with RCP8.5 \nscenario emissions (1765-1850)')
		ax1[0,0].fill_between(np.arange(1765,int(gcp_ems[-1,0])+1), CO2_ems_lower, CO2_ems_upper, color='orangered', alpha=0.2)
		ax1[0,0].set_xlim(1860,2020)

		# ax1[0,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), C_rcp[0,0,:], color='red')
		# ax1[0,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), C[0,0,:], color='green')
		ax1[0,1].fill_between(np.arange(1765,int(gcp_ems[-1,0])+1), C_fit_upper[0,0,:], C_fit_lower[0,0,:], color='orangered', alpha=0.2)
		ax1[0,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), C_fit[0,0,:], color='orange', linestyle='-')
		# ax1[0,1].plot(noaa_co2_annual_gl[:,0], noaa_co2_annual_gl[:,1], color='blue', linestyle='--')
		ax1[0,1].scatter(noaa_co2_annual_gl[:,0], noaa_co2_annual_gl[:,1], s=24, marker='o', color='black', alpha=0.5)
		ax1[0,1].set_xlabel('Year')
		ax1[0,1].set_ylabel('CO$_2$ concentration (ppmv)')
		ax1[0,1].set_xlim(1860,2020)

		# ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), RF_rcp[0,:], color='red')
		# ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), RF[0,:], color='green')
		ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), RF_fit[0,:], color='red', linestyle='-', label='total')
		ax1[1,0].fill_between(np.arange(1765,int(gcp_ems[-1,0])+1), RF_fit_lower[0,:] - other_rf, RF_fit_upper[0,:] - other_rf, color='orangered', alpha=0.2)
		ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), RF_fit[0,:] - other_rf, color='orange', linestyle='-', label='CO$_2$')
		ax1[1,0].plot(np.arange(1765,int(gcp_ems[-1,0])+1), other_rf, color='green', linestyle='-', label='non-CO$_2$')
		ax1[1,0].legend(loc='best', edgecolor='white', framealpha=0.)
		ax1[1,0].set_xlabel('Year')
		ax1[1,0].set_ylabel('Radiative forcing (Wm$^{-2}$)')
		ax1[1,0].set_xlim(1860,2020)

		ax1[1,0].plot(np.arange(1765,2018), other_rf_magicc[:2018-1765] + 5.78*np.log(C_fit[0,0,:]/278.), color='purple')
		T_magicc = multiscen_oxfair(emissions=CO2_ems, other_rf = other_rf_magicc[:2018-1765],emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),a = np.array([[0.2173,0.2240,0.2824,0.2763]]),tau = np.array([[1000000,394.4,36.54,4.304]]),r = np.array([[r0_fit[0],r0_fit[1],r0_fit[2],0.0]]),PI_C = np.array([278.0]),iirf100_max = 97.0,f = np.array([[5.35,0.,0.]]),tcr = tcr_val,ecs = ecs_val, d = np.array([239.0,4.1]),q = np.array([0.33,0.41]),F_2x = 3.74,multigas=False, multiscen=False)[2][0,:]
		ax1[1,1].plot(np.arange(1765,2018), T_magicc - np.mean(T_magicc[1850-1765:1901-1765]), color='purple', linestyle='-')

		# ax1[1,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), T_rcp[0,:] - np.mean(T_rcp[0,1850-1765:1901-1765]), color='red')
		# ax1[1,1].plot(np.arange(1765,int(gcp_ems[-1,0])+1), T[0,:] - np.mean(T[0,1850-1765:1901-1765]), color='green')
		ax1[1,1].set_xlabel('Year')
		ax1[1,1].set_ylabel('Temperature relative to 1850-1900 (K)')
		ax1[1,1].set_xlim(1860,2020)
		ax1[1,1].axhline(y = 1.0, linestyle=':', color='black', alpha=0.4)
		ax1[1,1].axhline(y = 0, linestyle='-', color='black', alpha=0.4)

		ax1[1,1].scatter(np.arange(1850,2018), annual_temp, marker='o', s=12, color='black', alpha=0.5, label='SR1.5 observations\n4-dataset-mean')
		ax1[1,1].legend(loc='best', edgecolor='white', framealpha=0.)


		T_other = multiscen_oxfair(emissions=np.zeros_like(other_rf), other_rf = other_rf,emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),a = np.array([[0.2173,0.2240,0.2824,0.2763]]),tau = np.array([[1000000,394.4,36.54,4.304]]),r = np.array([[r0_fit[0],r0_fit[1],r0_fit[2],0.0]]),PI_C = np.array([278.0]),iirf100_max = 97.0,f = np.array([[5.78,0.,0.]]),tcr = tcr_val,ecs = ecs_val, d = np.array([239.0,4.1]),q = np.array([0.33,0.41]),F_2x = 3.74,multigas=False, multiscen=False)[2][0,:]
		ax1[1,1].plot(np.arange(1765,2018), T_other - np.mean(T_other[1850-1765:1901-1765]), color='green', linestyle='-')

		scaled_rf = (piers_updated_rf[1765-1750:,14]-piers_updated_rf[1765-1750:,1])
		scaled_rf = np.append(scaled_rf, scaled_rf[-1])
		T_scaled = multiscen_oxfair(emissions=CO2_ems, other_rf = scaled_rf,emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),a = np.array([[0.2173,0.2240,0.2824,0.2763]]),tau = np.array([[1000000,394.4,36.54,4.304]]),r = np.array([[r0_fit[0],r0_fit[1],r0_fit[2],0.0]]),PI_C = np.array([278.0]),iirf100_max = 97.0,f = np.array([[5.78,0.,0.]]),tcr = tcr_val,ecs = ecs_val,d = np.array([239.0,4.1]),q = np.array([0.33,0.41]),F_2x = 3.74,multigas=False, multiscen=False)[2][0,:]
		ax1[1,1].plot(np.arange(1765,2018), T_scaled - np.mean(T_scaled[1850-1765:1901-1765]), color='orange', linestyle='-')

		T_resp = multiscen_oxfair(emissions=CO2_ems, other_rf = other_rf,emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),a = np.array([[0.2173,0.2240,0.2824,0.2763]]),tau = np.array([[1000000,394.4,36.54,4.304]]),r = np.array([[r0_fit[0],r0_fit[1],r0_fit[2],0.0]]),PI_C = np.array([278.0]),iirf100_max = 97.0,f = np.array([[5.78,0.,0.]]),tcr = tcr_val,ecs = ecs_val,d = np.array([239.0,4.1]),q = np.array([0.33,0.41]),F_2x = 3.74,multigas=False, multiscen=False)[2][0,:]
		ax1[1,1].plot(np.arange(1765,2018), T_resp - np.mean(T_resp[1850-1765:1901-1765]), color='red', linestyle='-')

		T_scaled_lower = multiscen_oxfair(emissions=CO2_ems_lower, other_rf = scaled_rf,emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),a = np.array([[0.2173,0.2240,0.2824,0.2763]]),tau = np.array([[1000000,394.4,36.54,4.304]]),r = np.array([[r0_fit_lower[0],r0_fit_lower[1],r0_fit_lower[2],0.0]]),PI_C = np.array([278.0]),iirf100_max = 97.0,f = np.array([[5.78,0.,0.]]),tcr = tcr_val,ecs = ecs_val,d = np.array([239.0,4.1]),q = np.array([0.33,0.41]),F_2x = 3.74,multigas=False, multiscen=False)[2][0,:]
		T_scaled_upper = multiscen_oxfair(emissions=CO2_ems_upper, other_rf = scaled_rf,emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.]) / 28.97),a = np.array([[0.2173,0.2240,0.2824,0.2763]]),tau = np.array([[1000000,394.4,36.54,4.304]]),r = np.array([[r0_fit_upper[0],r0_fit_upper[1],r0_fit_upper[2],0.0]]),PI_C = np.array([278.0]),iirf100_max = 97.0,f = np.array([[5.78,0.,0.]]),tcr = tcr_val,ecs = ecs_val,d = np.array([239.0,4.1]),q = np.array([0.33,0.41]),F_2x = 3.74,multigas=False, multiscen=False)[2][0,:]
		ax1[1,1].fill_between(np.arange(1765,2018), T_scaled_lower - np.mean(T_scaled_lower[1850-1765:1901-1765]), T_scaled_upper - np.mean(T_scaled_upper[1850-1765:1901-1765]), color='orangered', alpha=0.2)

		ax1[0,1].text(1870,410,'FIT -> r0 : %2.2f, rC : %1.4f, rT : %1.3f' % (r0_fit[0], r0_fit[1], r0_fit[2]))
		if x_val.success == True:
			ax1[0,1].text(1870,400,'Fit terminated successfully')
		else:
			ax1[0,1].text(1870,400,'Fit failed')

		# fig1.savefig('co2_fit.pdf', dpi=300)

		plt.show()



