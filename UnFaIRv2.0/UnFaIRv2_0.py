# UnFaIRv2.0

# structure

# 1 - Input format: Scenarios sets
	# input and output should be managed by dataframes 
	# use multiindexing in funciton to give user dataframe format to input emissions timeseries. 
	#     input['scen_name']['year']['gas_name'] = value

# 2 - input for parameter sets
    # use dataframe, where outer index is the parameter set number, inner index is the parameter name and column is the gas

# 3 - make into numpy array with correct dimensions

# 4 - compute output with standard functions

# 5 - format output in nice way, connecting emissions, C, RF, T and parameter set used in nice way 

import numpy as np
import pandas as pd

def return_empty_emissions(start_year=1765, end_year=2500, scen_names=[0]):

	df = pd.DataFrame({'CO2':np.zeros(end_year+1-start_year),'CH4':np.zeros(end_year+1-start_year),'N2O':np.zeros(end_year+1-start_year)}, index=np.arange(start_year,end_year+1))

	return pd.concat([df]*len(scen_names), keys=scen_names, axis=1)

def input_to_numpy(input_df):

	return input_df.values.reshape(input_df.columns.levels[0].size, input_df.columns.levels[1].size, input_df.index.size)


input_emissions = return_empty_emissions(1800,1820,np.arange(1,11))

input_numpy = input_to_numpy(input_emissions)

print(input_emissions)
print(input_numpy.shape)