import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

# returns thermal parameter dataframe of size param_num with varied params given by dictionary 'nonGas_params'
def thermal_param_df_creator(param_num, nonGas_params=None):

	standard_params = {'TCR':[1.75],'ECS':[2.6],'d_1':[4.1],'d_2':[239.0],'F_2x':[3.74]}
	parameters_df = pd.DataFrame(standard_params)
	parameters_df = parameters_df.reindex(range(0,param_num))

	if nonGas_params != None:
		for key in nonGas_params.keys():
			assert param_num == len(nonGas_params[key]), 'length of parameter %s input isn\'t same as param_num: %d, %d' % (key,len(nonGas_params[key]),param_num)
			assert key in parameters_df.keys(), 'param key "%s" isn\'t found in standard parameter name set' % (key)
			parameters_df[key] = nonGas_params[key]

	# fill NaNs with standard parameter values
	for key in standard_params.keys():
		parameters_df[key].fillna(standard_params[key][0], inplace=True)

	return parameters_df

# returns inputted parameter dataframe with columns for parameters of single gas
def add_new_gas_to_params(full_params_df, single_gas_params):

	assert full_params_df.index.values.shape == single_gas_params.index.values.shape, 'Index of two input dataframes are different'
	joined_df = pd.concat([full_params_df, single_gas_params], axis=1, sort=False)

	return joined_df

# returns standard parameters for gas with name gas_name ('CO2','N2O','CH4')
def standard_gas_param_df_creator(param_num, single_gas_name):

	if single_gas_name == 'CO2':
		standard_gas_params = {'f0_CO2':[3.74/np.log(2.0)],'f1_CO2':[0.0],'f2_CO2':[0.0],'iirf100_max_CO2':[97.0],'emis2conc_CO2':[28.97/(5.148*12.0)],'PI_C_CO2':[278.0],'r0_CO2':[32.40],'rC_CO2':[0.019],'rT_CO2':[4.165],'rA_CO2':[0.0],'a0_CO2':[0.2173],'a1_CO2':[0.2240],'a2_CO2':[0.2824],'a3_CO2':[0.2763],'tau0_CO2':[1e6],'tau1_CO2':[394.4],'tau2_CO2':[36.54],'tau3_CO2':[4.304]}
	elif single_gas_name == 'CH4':
		standard_gas_params = {'f0_CH4':[0.0],'f1_CH4':[0.0],'f2_CH4':[0.036],'iirf100_max_CH4':[97.0],'emis2conc_CH4':[28.97/(5.148*16.0)],'PI_C_CH4':[722.0],'r0_CH4':[9.05942806e+00],'rC_CH4':[-1.03745809e-07],'rT_CH4':[-1.85711888e-01],'rA_CH4':[1.45117387e-04],'a0_CH4':[1.0],'a1_CH4':[0.0],'a2_CH4':[0.0],'a3_CH4':[0.0],'tau0_CH4':[9.0],'tau1_CH4':[1.0],'tau2_CH4':[1.0],'tau3_CH4':[1.0]}
	elif single_gas_name == 'N2O':
		standard_gas_params = {'f0_N2O':[0.0],'f1_N2O':[0.0],'f2_N2O':[0.12],'iirf100_max_N2O':[97.0],'emis2conc_N2O':[28.97/(5.148*28.0)],'PI_C_N2O':[273.0],'r0_N2O':[4.97443512e+01],'rC_N2O':[5.87120814e-04],'rT_N2O':[-2.02130466e+00],'rA_N2O':[2.07719812e-02],'a0_N2O':[1.0],'a1_N2O':[0.0],'a2_N2O':[0.0],'a3_N2O':[0.0],'tau0_N2O':[121.0],'tau1_N2O':[1.0],'tau2_N2O':[1.0],'tau3_N2O':[1.0]}
	else:
		return 

	parameters_df = pd.DataFrame(standard_gas_params)
	parameters_df = parameters_df.reindex(range(0,param_num))
	# fill NaNs with standard parameter values
	for key in parameters_df.keys():
		parameters_df[key].fillna(parameters_df[key].loc[0], inplace=True)

	return parameters_df




param_num = 1000 # number of parameter sets 

TCR_vals = np.random.normal(loc=1.75,scale=0.3,size=param_num)
ECS_vals = np.random.normal(loc=2.6,scale=0.5,size=param_num)

params_to_vary = {'TCR':TCR_vals,'ECS':ECS_vals} # values to input into parameter dataframe



thermal_params_df = thermal_param_df_creator(param_num, nonGas_params=params_to_vary)

co2_standard_params = standard_gas_param_df_creator(param_num, 'CO2')
ch4_standard_params = standard_gas_param_df_creator(param_num, 'CH4')
n2o_standard_params = standard_gas_param_df_creator(param_num, 'N2O')

return_df = add_new_gas_to_params(thermal_params_df, co2_standard_params)
return_df = add_new_gas_to_params(return_df, ch4_standard_params)
return_df = add_new_gas_to_params(return_df, n2o_standard_params)

fig1 = return_df.TCR.hist(bins = 20)
plt.figure()
fig2 = return_df.ECS.hist(bins = 20)
plt.show()












