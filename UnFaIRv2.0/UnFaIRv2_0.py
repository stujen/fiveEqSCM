# UnFaIRv2.0 - Stuart Jenkins and Nick Leach

# structure

# 1 - Input format: Scenarios sets
	# input and output should be managed by dataframes
	# use multiindexing in funciton to give user dataframe format to input emissions timeseries.
	#     input['scen_name']['year']['gas_name'] = value

# 2 - input for parameter sets
    # use dataframe, where outer index is the parameter set number, inner index is the parameter name and column is the gas
	# potentially separate gas cycle parameters and thermal parameters since the model dimensions are different

# 3 - make into numpy array with correct dimensions
	# Dimensions : [scenario, gas params, thermal params, gas, time/gas pools]

# 4 - compute output with standard functions

# 5 - format output in nice way, connecting emissions, C, RF, T and parameter set used in nice way
	# Add help kwargs to assist users get the right input shape?

import numpy as np
import pandas as pd
import scipy as sp

def return_empty_emissions(start_year=1765, end_year=2500, scen_names=[0]):

	# returns an emissions dataframe of the correct format for use in UnFaIR with the given scenario names

	df = pd.DataFrame({'CO2':np.zeros(end_year+1-start_year),'CH4':np.zeros(end_year+1-start_year),'N2O':np.zeros(end_year+1-start_year)}, index=np.arange(start_year,end_year+1))

	df = pd.concat([df]*len(scen_names), keys=scen_names, axis=1)

	df.index = df.index.rename('Year')

	df.columns = df.columns.rename(['Scenario','Gas'])

	return df

def return_empty_forcing(start_year=1765, end_year=2500, scen_names=[0]):

	df = pd.DataFrame({'forcing':np.zeros(end_year+1-start_year)}, index=np.arange(start_year,end_year+1))

	df = pd.concat([df]*len(scen_names), keys=scen_names, axis=1)

	df.index = df.index.rename('Year')

	df.columns = df.columns.rename(['Scenario','Variable'])

	return df

	return

def input_to_numpy(input_df):

	# converts the dataframe input into a numpy array for calculation, dimension order = [name, gas, time/parameter]

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

	gas_cycle_parameters.index = gas_cycle_parameters.index.rename('param_name')

	gas_cycle_parameters.columns = gas_cycle_parameters.columns.rename(['Gas_cycle_set','Gas'])

	return gas_cycle_parameters.apply(pd.to_numeric)

def default_gas_forcing_param_uncertainty():

	# returns a dataframe of default parameters in the format UnFaIR requires (pd.concat -> additional sets)

	gas_parameter_list = ['a1','a2','a3','a4','tau1','tau2','tau3','tau4','r0','rC','rT','rA','PI_conc','emis2conc','f1','f2','f3']

	gas_parameter_uncertainty = pd.DataFrame(columns=['CO2','CH4','N2O'],index=gas_parameter_list)

	gas_parameter_uncertainty.loc['a1':'a4'] = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0]]).T
	gas_parameter_uncertainty.loc['tau1':'tau4'] = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0]]).T
	gas_parameter_uncertainty.loc['r0':'rA'] = np.array([[0.13,0.13,0.13,0],[0.10,0,0.14,0.19],[0.18,0,0,0.20]]).T
	gas_parameter_uncertainty.loc['PI_conc'] = np.array([0,0,0])
	gas_parameter_uncertainty.loc['emis2conc'] = np.array([0,0,0])
	gas_parameter_uncertainty.loc['f1':'f3'] = np.array([[0,0,0],[0,0,0],[0,0,0]]).T

	gas_parameter_uncertainty = pd.concat([gas_parameter_uncertainty], keys = ['normal'], axis = 1)

	gas_parameter_uncertainty.index = gas_parameter_uncertainty.index.rename('param_name')

	gas_parameter_uncertainty.columns = gas_parameter_uncertainty.columns.rename(['Distribution','Gas'])

	return gas_parameter_uncertainty.apply(pd.to_numeric)

def default_thermal_params():

	# returns a dataframe of default parameters in the format UnFaIR requires (pd.concat -> additional sets)

	thermal_parameter_list = ['d','tcr_ecs']

	thermal_parameters = pd.DataFrame(columns=[1,2],index=thermal_parameter_list)
	thermal_parameters.loc['d'] = np.array([239.0,4.1])
	thermal_parameters.loc['tcr_ecs'] = np.array([1.6,2.75])

	thermal_parameters = pd.concat([thermal_parameters], keys = ['default'], axis = 1)

	thermal_parameters.index = thermal_parameters.index.rename('param_name')

	thermal_parameters.columns = thermal_parameters.columns.rename(['Thermal_param_set','Box'])

	return thermal_parameters.apply(pd.to_numeric)

def default_thermal_param_uncertainty():

	# returns a dataframe of default parameters in the format UnFaIR requires (pd.concat -> additional sets)

	thermal_parameter_list = ['d','tcr_ecs']

	thermal_parameter_uncertainty = pd.DataFrame(columns=[1,2],index=thermal_parameter_list)
	thermal_parameter_uncertainty = pd.concat([thermal_parameter_uncertainty]*2,keys=[5,95],axis=1)
	thermal_parameter_uncertainty.loc['d'] = [239,1.6,239,8.4]
	thermal_parameter_uncertainty.loc['tcr_ecs'] = [1,1.6,2.5,4.5]

	thermal_parameter_uncertainty.index = thermal_parameter_uncertainty.index.rename('param_name')

	thermal_parameter_uncertainty.columns = thermal_parameter_uncertainty.columns.rename(['Percentile','Box'])

	return thermal_parameter_uncertainty.apply(pd.to_numeric)

def draw_monte_carlo_param_set(N , input_parameters , input_uncertainties , type = 'normal'):

	# function that takes a single set of parameter medians and corresponding % uncertainty dataframe, and creates a
	# new dataframe with N samples of the parameter set.

	if type == 'normal':

		param_set = [input_parameters[input_parameters.columns.levels[0][0]]]

		for i in np.arange(N):

			param_set += [param_set[0] * np.random.normal(np.ones(input_parameters.shape),input_uncertainties)]

		param_set = pd.concat(param_set, keys = ['median']+[x + type for x in [str(i) for i in np.arange(N-1)]], axis = 1)

		return param_set

	if type == 'lognormal':

		param_set = input_parameters[input_parameters.columns.levels[0][0]]

		loc = ((param_set**2 - input_uncertainties[5]*input_uncertainties[95]) / (input_uncertainties[5]+input_uncertainties[95]-2*param_set)).fillna(0)
		mu = np.log(param_set+loc)
		scale = ( np.log(input_uncertainties[95]+loc) - mu ) / 1.645

		# Constrain to be within +/- 3 sigma
		constrain_high = param_set.copy()
		constrain_low = param_set.copy()
		constrain_low.loc[:] = sp.stats.lognorm.ppf(0.003,scale.fillna(0),-loc.fillna(0),np.exp(mu.fillna(0)))
		constrain_high.loc[:] = sp.stats.lognorm.ppf(0.973,scale.fillna(0),-loc.fillna(0),np.exp(mu.fillna(0)))
		constrain_low,constrain_high = constrain_low.fillna(-10**10),constrain_high.fillna(10**10)

		param_set = [param_set]

		for i in np.arange(N):
		    while True:
		        new_param_set = np.random.lognormal(mu.fillna(0),scale.fillna(0))-loc.fillna(0)

		        if all(new_param_set<constrain_high) and all(new_param_set>constrain_low):

		            break

		    param_set += [new_param_set]

		param_set = pd.concat(param_set, keys = ['median']+[x + 'lognorm' for x in [str(i) for i in np.arange(N-1)]], axis = 1)

		return param_set

def tcr_ecs_to_q(input_parameters=True , F_2x=3.74 , help=False):

	# converts a tcr / ecs / d dataframe into a d / q dataframe for use in UnFaIRv2

	if help:
		tcr_ecs_test = default_thermal_params()
		tcr_ecs_test = pd.concat([tcr_ecs_test['default']]*2,keys=['default','1'],axis=1)
		tcr_ecs_test.loc['tcr_ecs'] = [1.6,2.75,1.4,2.4]
		tcr_ecs_test = tcr_ecs_test.loc[['d','tcr_ecs']]
		print('Example input format:')
		return tcr_ecs_test

	if type(input_parameters.columns) != pd.core.indexes.multi.MultiIndex:
		return 'input_parameters not in MultiIndex DataFrame. Set help=True for formatting of input.'
	else:
		output_params = input_parameters.copy()
		param_arr = input_to_numpy(input_parameters)
		k = 1.0 - (param_arr[:,:,0]/70.0)*(1.0 - np.exp(-70.0/param_arr[:,:,0]))
		output_params.loc['q'] = ( ( param_arr[:,0,1][:,np.newaxis] - param_arr[:,1,1][:,np.newaxis] * np.roll(k,shift=1) )/( F_2x * ( k - np.roll(k,shift=1) ) ) ) .flatten()

		return output_params.loc[['d','q']]

def calculate_alpha(G,G_A,T,r,g0,g1,iirf100_max = 97.0):

	iirf100_val = r[...,0] + r[...,1] * (G-G_A) + r[...,2] * T + r[...,3] * G_A

	iirf100_val = np.abs(iirf100_val)

	iirf100_val = (iirf100_val>iirf100_max) * iirf100_max + iirf100_val * (iirf100_val<iirf100_max)

	alpha_val = g0 * np.sinh(iirf100_val / g1)

	return alpha_val

def step_concentration(R,E,alpha,a,tau,PI_conc,emis2conc):

	R = E * emis2conc[...,np.newaxis] * a * alpha * tau * ( 1. - np.exp( -1./(alpha*tau) ) ) + R * np.exp( -1./(alpha * tau) )

	C = PI_conc + np.sum(R,axis=-1)

	G_A = (C - PI_conc) / emis2conc

	return C,R,G_A

def step_forcing(C,PI_conc,f):

	RF = f[...,0] * np.log( C/PI_conc ) + f[...,1] * ( C - PI_conc ) + f[...,2] * ( np.sqrt(C) - np.sqrt(PI_conc) )

	return RF

def step_temperature(S,F,q,d):

	S = q * F * ( 1 - np.exp(-1/d) ) + S * np.exp(-1/d)

	T = np.sum(S,axis=-1)

	return S,T

# Run modes:
	# Full forward
	# Concentration driven
	# Set temperature response (gas cycle mode)
# Checks on:
	# timeseries length
	# parameter formatting
	# parameter size / shape
	# same number of scenarios in emissions and forcing
	#

def run_UnFaIR( emissions_in , \
				forcing_in = 0.0 , \
				gas_parameters = default_gas_forcing_params() , \
				thermal_parameters = default_thermal_params() , \
				show_run_info = True):

	# Determine the number of scenario runs , parameter sets , gases , integration period

	dim_scenario = emissions_in.columns.levels[0].size
	scen_names = list(emissions_in.columns.levels[0])
	dim_gas_param = gas_parameters.columns.levels[0].size
	gas_set_names = list(gas_parameters.columns.levels[0])
	dim_thermal_param = thermal_parameters.columns.levels[0].size
	thermal_set_names = list(thermal_parameters.columns.levels[0])
	n_gas = emissions_in.columns.levels[1].size
	n_year = emissions_in.index.size

	if show_run_info:
		print('Integrating ' + str(dim_scenario) + ' scenarios, ' + \
			   str(dim_gas_param) + ' gas cycle parameter sets, ' + \
			   str(dim_thermal_param) + ' thermal response parameter sets, over ' + \
			   str(list(emissions_in.columns.levels[1])) + ', between ' + \
			   str(emissions_in.index[0]) + ' and ' + str(emissions_in.index[-1]) + '...')

	# Slice the parameter sets into numpy arrays of the right shape
	# Dimensions : [scenario, gas params, thermal params, gas, time/gas pools]

	a = input_to_numpy(gas_parameters.loc['a1':'a4'])[np.newaxis,:,np.newaxis,...]
	tau = input_to_numpy(gas_parameters.loc['tau1':'tau4'])[np.newaxis,:,np.newaxis,...]
	r = input_to_numpy(gas_parameters.loc['r0':'rA'])[np.newaxis,:,np.newaxis,...]
	PI_conc = gas_parameters.loc['PI_conc'].values.reshape(gas_parameters.loc['PI_conc'].index.levels[0].size,gas_parameters.loc['PI_conc'].index.levels[1].size)[np.newaxis,:,np.newaxis,...]
	emis2conc = gas_parameters.loc['emis2conc'].values.reshape(gas_parameters.loc['emis2conc'].index.levels[0].size,gas_parameters.loc['emis2conc'].index.levels[1].size)[np.newaxis,:,np.newaxis,...]

	f = input_to_numpy(gas_parameters.loc['f1':'f3'])[np.newaxis,:,np.newaxis,...]

	thermal_parameters = tcr_ecs_to_q(thermal_parameters)

	d = thermal_parameters.loc['d'].values.reshape(thermal_parameters.loc['d'].index.levels[0].size,thermal_parameters.loc['d'].index.levels[1].size)[np.newaxis,np.newaxis,...]
	q = thermal_parameters.loc['q'].values.reshape(thermal_parameters.loc['q'].index.levels[0].size,thermal_parameters.loc['q'].index.levels[1].size)[np.newaxis,np.newaxis,...]

	# Reformat inputs into the right shape

	emissions = input_to_numpy(emissions_in)[:,np.newaxis,np.newaxis,...]

	if type(forcing_in)==float or type(forcing_in)==int:
		ext_forcing = forcing_in + input_to_numpy(return_empty_forcing(emissions_in.index[0],emissions_in.index[-1],emissions_in.columns.levels[0]))[:,np.newaxis,np.newaxis,...]
	elif type(forcing_in)==pd.core.frame.DataFrame:
		if type(forcing_in.columns)==pd.core.indexes.multi.MultiIndex:
			if forcing_in.index.equals(emissions_in.index):
				ext_forcing = input_to_numpy(forcing_in)[:,np.newaxis,np.newaxis,...]
			else:
				return "forcing timeseries length different to emission timeseries"
		else:
			return "forcing DataFrame not MultiIndex, use pd.concat([df],keys=['scenario1'...],axis=1)"
	else:
		return "forcing not pandas DataFrame, use return_empty_forcing to check correct inpurt formatting"


	# Create appropriate shape variable arrays

	G = np.cumsum(emissions,axis=-1)
	C = np.zeros((dim_scenario,dim_gas_param,dim_thermal_param,n_gas,n_year))
	RF = np.zeros((dim_scenario,dim_gas_param,dim_thermal_param,n_gas,n_year))
	T = np.zeros((dim_scenario,dim_gas_param,dim_thermal_param,n_year))
	alpha = np.zeros((dim_scenario,dim_gas_param,dim_thermal_param,n_gas,n_year))

	# Initialize the first timestep

	g1 = np.sum( a * tau * ( 1. - ( 1. + 100/tau ) * np.exp(-100/tau) ), axis=-1 )
	g0 = ( np.sinh( np.sum( a * tau * ( 1. - np.exp(-100/tau) ) , axis=-1) / g1 ) )**(-1.)

	alpha[...,0] = calculate_alpha(G=np.zeros(C[...,0].shape),G_A=np.zeros(C[...,0].shape),T=T[...,0,np.newaxis],r=r,g0=g0,g1=g1)
	C[...,0],R,G_A = step_concentration(R = np.zeros(a.shape),alpha=alpha[...,0,np.newaxis],E=emissions[...,0,np.newaxis],\
										a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)
	RF[...,0] = step_forcing(C=C[...,0],PI_conc=PI_conc,f=f)
	S,T[...,0] = step_temperature(S=np.zeros(d.shape),F=np.sum(RF[...,0],axis=-1)[...,np.newaxis]+ext_forcing[...,0],q=q,d=d)

	# Step over remaining timesteps

	for t in np.arange(1,emissions.shape[-1]):

		alpha[...,t] = calculate_alpha(G=G[...,t-1],G_A=G_A,T=T[...,t-1,np.newaxis],r=r,g0=g0,g1=g1)
		C[...,t],R,G_A = step_concentration(R = R,alpha=alpha[...,t,np.newaxis],E=emissions[...,t,np.newaxis],\
											a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)
		RF[...,t] = step_forcing(C=C[...,t],PI_conc=PI_conc,f=f)
		S,T[...,t] = step_temperature(S=S,F=np.sum(RF[...,t],axis=-1)[...,np.newaxis]+ext_forcing[...,t],q=q,d=d)

	RF = np.concatenate((RF,np.sum(RF,axis=-2)[...,np.newaxis,:]),axis=-2)

	out_dict = { \
				'C':pd.DataFrame(C.T.swapaxes(1,-1).swapaxes(2,-2).reshape(n_year,n_gas*dim_scenario*dim_gas_param*dim_thermal_param),index = emissions_in.index,columns=pd.MultiIndex.from_product([scen_names,gas_set_names,thermal_set_names,['CO2','CH4','N2O']],names=['Scenario','Gas cycle set','Thermal set','Gas name'])), \
				'RF':pd.DataFrame(RF.T.swapaxes(1,-1).swapaxes(2,-2).reshape(n_year,(n_gas+1)*dim_scenario*dim_gas_param*dim_thermal_param),index = emissions_in.index,columns=pd.MultiIndex.from_product([scen_names,gas_set_names,thermal_set_names,['CO2','CH4','N2O','Total']],names=['Scenario','Gas cycle set','Thermal set','Gas name'])), \
				'T':pd.DataFrame(T.T.swapaxes(1,-1).reshape(n_year,dim_scenario*dim_gas_param*dim_thermal_param),index = emissions_in.index,columns=pd.MultiIndex.from_product([scen_names,gas_set_names,thermal_set_names],names=['Scenario','Gas cycle set','Thermal set'])), \
				'alpha':pd.DataFrame(alpha.T.swapaxes(1,-1).swapaxes(2,-2).reshape(n_year,n_gas*dim_scenario*dim_gas_param*dim_thermal_param),index = emissions_in.index,columns=pd.MultiIndex.from_product([scen_names,gas_set_names,thermal_set_names,['CO2','CH4','N2O']],names=['Scenario','Gas cycle set','Thermal set','Gas name'])), \
			   }

	for axis in out_dict.keys():
		out_dict[axis].index = out_dict[axis].index.rename('Year')

	return out_dict
