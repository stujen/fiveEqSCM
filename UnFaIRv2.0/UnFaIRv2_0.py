# UnFaIRv2.0

# structure

# 1 - Input format: Scenarios sets
	# input and output should be managed by dataframes
	# use multiindexing in funciton to give user dataframe format to input emissions timeseries.
	#     input['scen_name']['year']['gas_name'] = value

# 2 - input for parameter sets
    # use dataframe, where outer index is the parameter set number, inner index is the parameter name and column is the gas
	# potentially separate gas cycle parameters and thermal parameters since the model dimensions are different

# 3 - make into numpy array with correct dimensions
	# Dimensions : [scenario, gas params, thermal params, gas, time]

# 4 - compute output with standard functions

# 5 - format output in nice way, connecting emissions, C, RF, T and parameter set used in nice way

import numpy as np
import pandas as pd

def return_empty_emissions(start_year=1765, end_year=2500, scen_names=[0]):

	# returns an emissions dataframe of the correct format for use in UnFaIR with the given scenario names

	df = pd.DataFrame({'CO2':np.zeros(end_year+1-start_year),'CH4':np.zeros(end_year+1-start_year),'N2O':np.zeros(end_year+1-start_year)}, index=np.arange(start_year,end_year+1))

	return pd.concat([df]*len(scen_names), keys=scen_names, axis=1)

def input_to_numpy(input_df):

	# converts the dataframe input into a numpy array for calculation, dimension order = [scenario, gas, time]

	return input_df.values.T.reshape(input_df.columns.levels[0].size, input_df.columns.levels[1].size, input_df.index.size)

def default_gas_forcing_params():

	# returns a dataframe of default parameters in the format UnFaIR requires (pd.concat -> additional sets)

	gas_parameter_list = ['a1','a2','a3','a4','tau1','tau2','tau3','tau4','r0','rC','rT','rA','PI_conc','emis2conc','f1','f2','f3']

	gas_cycle_parameters = pd.DataFrame(columns=['CO2','CH4','N2O'],index=gas_parameter_list)

	gas_cycle_parameters.loc['a1':'a4'] = np.array([[0.2173,0.2240,0.2824,0.2763],[1,0,0,0],[1,0,0,0]]).T
	gas_cycle_parameters.loc['tau1':'tau4'] = np.array([[1000000,394.4,36.54,4.304],[9.15,1,1,1],[116.,1,1,1]]).T
	gas_cycle_parameters.loc['r0':'rA'] = np.array([[29.5,0.018,3.946,0],[9.15,0,-0.241,0.000263],[61.0,0,0,-0.000648]]).T
	gas_cycle_parameters.loc['PI_conc'] = np.array([278.0,700.0,276.0])
	gas_cycle_parameters.loc['emis2conc'] = 1/(5.148*10**18/1e18*np.array([12.,16.,28.])/28.97)
	gas_cycle_parameters.loc['f1':'f3'] = np.array([[5.78188211,0,0],[0,0,0.03895942],[0,0,0.11082109]]).T

	gas_cycle_parameters = pd.concat([gas_cycle_parameters], keys = ['default'], axis = 1)

	return gas_cycle_parameters.apply(pd.to_numeric)

def default_thermal_params():

	# returns a dataframe of default parameters in the format UnFaIR requires (pd.concat -> additional sets)

	thermal_parameter_list = ['d','q']

	thermal_parameters = pd.DataFrame(columns=[1,2],index=thermal_parameter_list)
	thermal_parameters.loc['d'] = np.array([239.0,4.1])
	thermal_parameters.loc['q'] = np.array([0.33,0.41])

	thermal_parameters = pd.concat([thermal_parameters], keys = ['default'], axis = 1)

	return thermal_parameters.apply(pd.to_numeric)



# input_emissions = return_empty_emissions(1800,1820,np.arange(1,11))
#
# input_numpy = input_to_numpy(input_emissions)
#
# print(input_emissions)
# print(input_numpy.shape)
