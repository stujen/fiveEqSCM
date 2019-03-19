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



