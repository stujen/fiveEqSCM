# A base version of the Universal Fractionally adjusted Impulse Response model (UnFaIR)

# Written by Stuart Jenkins before Nick Leach got hold of it...

# Documentation format:

# name (shape) description

## Dependencies ##

import numpy as np
import string
import math
import pandas as pd
import scipy as sp

def g_1(a,tau,h=100):

    # Calculates the g1 coefficient for the alpha function

    # a (3,4) gas cycle pool fractions
    # tau (3,4) gas cycle pool decay timescales
    # h (1) impulse response time horizon

    g1 = np.sum( a * tau * ( 1. - ( 1. + h/tau ) * np.exp(-h/tau) ), axis=1 )

    return g1

def g_0(a,tau,h=100):

    # Calculates the g0 coefficient for the alpha function

    # a (3,4) gas cycle pool fractions
    # tau (3,4) gas cycle pool decay timescales
    # h (1) impulse response time horizon

    g0 = ( np.sinh( np.sum( a * tau * ( 1. - np.exp(-h/tau) ) , axis=1) / g_1(a,tau,h) ) )**(-1.)

    return g0

def alpha_val(G,G_A,T,a,tau,r,h=100,iirf100_max = 97.0):

    # Computes alpha based on given parameters and variables for a single timestep

    # r (3,4) gas cycle parameters [r0 , rC , rT , rA]
    # T (1) current temperature
    # G (1) cumulative emissions
    # G_A (1) atmospheric gas burden in emission units
    # h (1) impulse response time horizon

    iirf100_val = r[...,0] + \
                  r[...,1] * (G-G_A) + \
                  r[...,2] * T + \
                  r[...,3] * G_A

    iirf100_val = np.abs(iirf100_val)

    iirf100_val = (iirf100_val>iirf100_max) * iirf100_max + iirf100_val * (iirf100_val<iirf100_max)

    alpha_val = g_0(a,tau,h) * np.sinh(iirf100_val / g_1(a,tau,h))

    return alpha_val

def step_conc(R,E,alpha,a,tau,PI_conc,emis2conc):

    # Computes concentrations from emissions for a single timestep

    # alpha (3) timescale adjustment factor
    # R (3,4) gas cycle pool concentrations
    # E (3) emissions
    # emis2conc (3) conversion factor from emission -> concentration units
    # a (3,4) gas cycle pool fractions
    # tau (3,4) gas cycle pool timescales
    # C (3) current concentration
    # PI_conc (3) pre industrial concentration
    # G_A (3) atmospheric gas burden expressed in emission units


    alpha = alpha[:,np.newaxis]

    R = E[:,np.newaxis] * emis2conc[:,np.newaxis] * a * alpha * tau * ( 1. - np.exp( -1./(alpha*tau) ) ) + R * np.exp( -1./(alpha * tau) )

    C = PI_conc + np.sum(R,axis=1)

    G_A = (C - PI_conc) / emis2conc

    return C,R,G_A

def step_forc(C,PI_conc,f):

    # Computes forcing from concentrations for a single timestep

    # C (3) current concentration
    # PI_conc (3) preindustrial concentration
    # F_ext (1) external forcing
    # f (3,3) gas forcing parameters

    RF = f[...,0] * np.log( C/PI_conc ) + \
        f[...,1] * ( C - PI_conc ) + \
        f[...,2] * ( np.sqrt(C) - np.sqrt(PI_conc) )

    return RF

def step_temp(S,F,q,d):

    # Computes temperature from total forcing for a single timestep

    # S (2) thermal boxes
    # F (1) current forcing
    # q (2) thermal box heat capacities
    # d (2) thermal box decay timescales

    S = q * F * ( 1 - np.exp(-1/d) ) + S * np.exp(-1/d)

    T = np.sum(S)

    return S,T

def default_gas_params():

    # Function that returns the default gas cycle parameters in a dataframe

    gas_param_list = ['a1','a2','a3','a4','tau1','tau2','tau3','tau4','r0','rC','rT','rA','PI_conc','emis2conc']
    gas_cycle_parameters = pd.DataFrame(columns=['CO2','CH4','N2O'],index=gas_param_list)
    gas_cycle_parameters.loc['a1':'a4'] = np.array([[0.2173,0.2240,0.2824,0.2763],[1,0,0.,0.],[1,0,0.,0.]]).T
    gas_cycle_parameters.loc['tau1':'tau4'] = np.array([[1000000,394.4,36.54,4.304],[9.15,1,1,1],[116.,1,1,1]]).T
    gas_cycle_parameters.loc['r0':'rA'] = np.array([[37.493303,0.01909,3.616153,0.0],\
                  [ 8.540000, 0, -0.360000,  0.000310],\
                  [ 67.231092,  0, 0,  -0.000906]]).T
    gas_cycle_parameters.loc['PI_conc'] = np.array([278.0,700.0,273.0])
    gas_cycle_parameters.loc['emis2conc'] = 1/(5.148*10**18 / 1e18 * np.array([12.,16.,28.]) / 28.97)

    return gas_cycle_parameters.apply(pd.to_numeric)

def default_forcing_params():

    # Function that returns the default forcing parameters in a dataframe

    forcing_param_list = ['f1','f2','f3']
    forcing_params = pd.DataFrame(columns=['CO2','CH4','N2O'],index=forcing_param_list)
    forcing_params.loc['f1':'f3'] = np.array([[3.74/np.log(2.),0.,0.],[0,0.,0.036],[0,0,0.12]]).T

    return forcing_params.apply(pd.to_numeric)

def default_thermal_params():

    # Function that returns the default thermal parameters in a dataframe

    thermal_param_list = ['d','q']
    thermal_params = pd.DataFrame(columns=[1,2],index=thermal_param_list)
    thermal_params.loc['d'] = np.array([239.0,4.1])
    thermal_params.loc['q'] = np.array([0.33,0.41])

    return thermal_params.apply(pd.to_numeric)

def tcrecs_to_q(tcr, \
                ecs, \
                d = np.array([239.0,4.1]), \
                F_2x = 3.74):

    k = 1.0 - (d/70.0)*(1.0 - np.exp(-70.0/d))

    if (tcr and ecs):

        q =  ( tcr - ecs * np.roll(k,shift=1) )/( F_2x * ( k - np.roll(k,shift=1) ) )

    return q

def UnFaIR(emissions_in, \
           F_ext = 0, \
           gas_params = default_gas_params(), \
           forcing_params = default_forcing_params(), \
           thermal_params = default_thermal_params() ):

    # The UnFaIR model... will not check input shapes are consistent

    # Slice relevant values from parameter dataframes:

    a = gas_params.loc['a1':'a4'].values.T
    tau = gas_params.loc['tau1':'tau4'].values.T
    r = gas_params.loc['r0':'rA'].values.T
    PI_conc = gas_params.loc['PI_conc'].values
    emis2conc = gas_params.loc['emis2conc'].values

    f = forcing_params.values.T

    d = thermal_params.loc['d'].values
    q = thermal_params.loc['q'].values

    # Reformat emissions_in into a numpy array:

    emissions = emissions_in.values.T

    # Set up other variable arrays:

    G = np.cumsum(emissions,axis=1)
    C = np.zeros(emissions.shape)
    RF = np.zeros(emissions.shape)
    T = np.zeros(emissions[0].shape)
    alpha = np.zeros(emissions.shape)

    if type(F_ext) == np.ndarray:
        if F_ext.size != emissions[0].size:
            return 'ERROR : External forcing length does not equal emissions input length'

    else:
        F_ext = np.full(emissions[0].size,float(F_ext))

    # Initialize the first timestep

    alpha[...,0] = alpha_val(G=0,G_A=0,T=0,tau=tau,a=a,r=r)
    C[...,0],R,G_A = step_conc(R = np.zeros((3,4)),alpha=alpha[...,0],E=emissions[...,0],a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)
    RF[...,0] = step_forc(C=C[...,0],PI_conc=PI_conc,f=f)
    S,T[0] = step_temp(S=np.zeros(2),F=np.sum(RF[...,0])+F_ext[0],q=q,d=d)

    # Step over the remaining times

    for t in np.arange(1,emissions[0].size):

        alpha[...,t] = alpha_val(G=G[...,t-1],G_A=G_A,T=T[t-1],tau=tau,a=a,r=r)
        C[...,t],R,G_A = step_conc(R = R,alpha=alpha[...,t],E=emissions[...,t],a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)
        RF[...,t] = step_forc(C=C[...,t],PI_conc=PI_conc,f=f)
        S,T[t] = step_temp(S=S,F=np.sum(RF[...,t])+F_ext[t],q=q,d=d)

    # Create the result output

    result = {'C': pd.DataFrame(C.T, index = emissions_in.index, columns = emissions_in.columns), \
              'RF': pd.DataFrame(np.concatenate((RF.T,F_ext[:,np.newaxis]),axis=1), index = emissions_in.index, columns = list(emissions_in.columns)+['F_ext']), \
              'T': pd.DataFrame(T.T, index = emissions_in.index, columns = ['Total']), \
              'alpha': pd.DataFrame(alpha.T, index = emissions_in.index, columns = emissions_in.columns), \
              'E':emissions_in, \
              'gas_params':gas_params, \
              'forcing_params':forcing_params, \
              'thermal_params':thermal_params}

    result['RF']['Total'] = result['RF'].sum(axis=1)

    return result

def fit_gas_cycles(emissions_in, \
                   temp_input, \
                   gas_params = default_gas_params()):


    a = gas_params.loc['a1':'a4'].values.T
    tau = gas_params.loc['tau1':'tau4'].values.T
    r = gas_params.loc['r0':'rA'].values.T
    PI_conc = gas_params.loc['PI_conc'].values
    emis2conc = gas_params.loc['emis2conc'].values

    emissions=emissions_in.values.T

    G = np.cumsum(emissions,axis=1)
    C = np.zeros(emissions.shape)
    alpha = np.zeros(emissions.shape)

    T = temp_input.values.flatten()

    alpha[...,0] = alpha_val(G=0,G_A=0,T=0,tau=tau,a=a,r=r)
    C[...,0],R,G_A = step_conc(R = np.zeros((3,4)),alpha=alpha[...,0],E=emissions[...,0],a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)

    for t in np.arange(1,emissions[0].size):

        alpha[...,t] = alpha_val(G=G[...,t-1],G_A=G_A,T=T[t-1],tau=tau,a=a,r=r)
        C[...,t],R,G_A = step_conc(R = R,alpha=alpha[...,t],E=emissions[...,t],a=a,tau=tau,PI_conc=PI_conc,emis2conc=emis2conc)

    results = pd.DataFrame(C.T, index = emissions_in.index, columns = emissions_in.columns)

    return results

def plot_all(outputs,names=['0','1','2','3']):
    from matplotlib import pyplot as plt
    import matplotlib

    color_dict = {'CO2':'r','CH4':'b','N2O':'g','Total':'k','F_ext':'brown'}

    spec_list = ['CO2','CH4','N2O','Total','F_ext']

    unit_dict = {'E': {'CO2':'GtC','CH4':'GtCH$_4$','N2O':'GtN'},\
                 'C': {'CO2':'ppm','CH4':'ppb','N2O':'ppb'},\
                 'alpha': {'CO2':'-','CH4':'yrs','N2O':'yrs'},\
                 'All':{'RF':'Wm$^{-2}$'} ,'Total':{'T':'Deg C'}}

    linestyles = ['-','--','-.',':']

    fig = plt.figure(figsize=(15,15))
    gs = matplotlib.gridspec.GridSpec(4,6)
    ax={}
    ax['CO2'] = {}
    ax['CH4'] = {}
    ax['N2O'] = {}
    ax['RF'] = {}
    ax['T'] = {}
    ax['CO2']['E'] = plt.subplot(gs[0,:2])
    ax['CH4']['E'] = plt.subplot(gs[0,2:4])
    ax['N2O']['E'] = plt.subplot(gs[0,4:])
    ax['CO2']['C'] = plt.subplot(gs[1,:2])
    ax['CH4']['C'] = plt.subplot(gs[1,2:4])
    ax['N2O']['C'] = plt.subplot(gs[1,4:])
    ax['CO2']['alpha'] = plt.subplot(gs[2,:2])
    ax['CH4']['alpha'] = plt.subplot(gs[2,2:4])
    ax['N2O']['alpha'] = plt.subplot(gs[2,4:])
    ax['RF']['All'] = plt.subplot(gs[3,:3])
    ax['T']['Total'] = plt.subplot(gs[3,3:])

    for i,output in enumerate(outputs):

        # Plot emissions:

        for s in spec_list[:3]:
            output['E'][s].plot(ax=ax[s]['E'],color = color_dict[s],linestyle=linestyles[i],label=names[i])

        # Plot concentrations:
        for s in spec_list[:3]:
            output['C'][s].plot(ax=ax[s]['C'],color = color_dict[s],linestyle=linestyles[i],label=names[i])

        # Plot alpha for CO2 / tau value for CH4/N2O:
        output['alpha']['CO2'].plot(ax=ax['CO2']['alpha'],color=color_dict['CO2'],linestyle=linestyles[i],label=names[i])
        for s in spec_list[1:3]:
            (output['alpha'][s]*output['gas_params'].loc['tau1'][s]).plot(ax=ax[s]['alpha'],color = color_dict[s],linestyle=linestyles[i],label=names[i])

        # Plot RFs:
        for s in spec_list:
            output['RF'][s].plot(ax=ax['RF']['All'],color=color_dict[s],label=s+'_'+names[i],linestyle=linestyles[i])

        #Plot Temp
        output['T']['Total'].plot(ax=ax['T']['Total'],color='k',label=names[i],linestyle=linestyles[i])

    #Label
    for key1 in ax.keys():
        for key2 in ax[key1].keys():
            ax[key1][key2].set_title(key1 + '_' + key2)
            ax[key1][key2].set_xlabel('Year')
            ax[key1][key2].set_ylabel(unit_dict[key2][key1])
            ax[key1][key2].legend()

    ax['CH4']['alpha'].set_title('lifetime')
    ax['N2O']['alpha'].set_title('lifetime')

    plt.tight_layout()


def Unstep_temp(Temp , thermal_params = default_thermal_params()):

    # Temp should be a panda series

    q = thermal_params.loc['q'].values
    d = thermal_params.loc['d'].values

    #Set up required arrays:

    T = Temp.values.flatten()
    A = q * ( 1 - np.exp(-1/d) )
    F = T.copy()

    # Initialize @ t=0

    F[0] = T[0] / np.sum(A)

    # Loop over timesteps:

    for t in np.arange(1,T.size):
        F[t] = (1/np.sum(A))*(T[t] - np.sum(A * F[:t][:,np.newaxis] * np.exp((np.arange(t)[:,np.newaxis] - t)/d)))

    return {'RF': pd.Series(data=F,index = Temp.index),\
            'T_in': Temp, \
            'thermal_params': thermal_params}

def Unstep_forc(RF , F_ext = 0 , forcing_params = default_forcing_params() , gas_params = default_gas_params()):

    # RF should be a panda series
    # F_ext should be a constant or array of the same length of RF

    if type(F_ext) == np.ndarray:
        if F_ext.size != RF.size:
            return 'ERROR : External forcing length does not equal emissions input length'

    f = forcing_params['CO2'].values
    PI_conc = gas_params.loc['PI_conc']['CO2']

    if ~f[1:].any():
        return {'CO2-fe_C': pd.Series(PI_conc * np.exp( ( RF.values - F_ext ) / f[0] ) , index = RF.index), \
                'RF_in': pd.Series(RF.values - F_ext, index=RF.index), \
                'forcing_params': forcing_params , 'gas_params' : gas_params}

    def num_solve(x):

        return step_forc(x, PI_conc, f) - (RF.values - F_ext)

    return {'CO2-fe_C': pd.Series(sp.optimize.root(num_solve, np.ones(RF.size), method='hybr')['x'], index=RF.index), \
            'RF_in': pd.Series(RF.values - F_ext, index=RF.index), \
            'forcing_params': forcing_params , 'gas_params' : gas_params}

def Unstep_concs(C , T , gas_params = default_gas_params()):

    # C & T must be panda series of the same length

    a = gas_params.loc['a1':'a4'].CO2.values[np.newaxis,:]
    tau = gas_params.loc['tau1':'tau4'].CO2.values[np.newaxis,:]
    r = gas_params.loc['r0':'rA'].CO2.values[np.newaxis,:]
    PI_conc = gas_params.loc['PI_conc'].CO2
    emis2conc = gas_params.loc['emis2conc'].CO2

    if not(T.index.equals(C.index)):
        return 'Temperature and concentration arrays are not the same size OR timeperiod!'

    T_vals = T.values
    C_vals = C.values

    alpha = np.zeros(T_vals.size)
    G = alpha.copy()
    G_A = alpha.copy()
    emissions = alpha.copy()

    G_A = (C_vals - PI_conc) / emis2conc

    alpha[0] = alpha_val(0,0,0,a,tau,r)

    for t in np.arange(1,T_vals.size):

        G = np.cumsum(emissions)

        alpha[t] = alpha_val(G[t-1],G_A[t-1],T_vals[t-1],a,tau,r)

        emissions[t] = ( (C_vals[t] - PI_conc)/emis2conc - \
                        np.sum(emissions[:t,np.newaxis]*a*tau*alpha[:t,np.newaxis]*(1-np.exp(-1/(alpha[:t,np.newaxis]*tau))) * \
                        np.cumprod(np.exp(-1 / (alpha[:t,np.newaxis]*tau))[::-1],axis=0)[::-1]) ) / \
                        (np.sum(a*tau*alpha[t]*(1-np.exp(-1/(alpha[t]*tau))))*emis2conc)

    return {'CO2-fe_E': pd.Series(data=emissions, index = T.index) ,\
            'alpha': pd.Series(data=alpha,index=T.index),\
            'C_in' : C, 'T_in' : T, \
            'gas_params': gas_params}

def UnUnFaIR(T, F_ext= 0 , gas_params = default_gas_params(), thermal_params= default_thermal_params(), forcing_params= default_forcing_params()):

    RF = Unstep_temp( T, thermal_params=thermal_params )['RF']
    C = Unstep_forc(RF, F_ext = F_ext, forcing_params=forcing_params, gas_params=gas_params)['CO2-fe_C']
    E = Unstep_concs(C, T, gas_params=gas_params)['CO2-fe_E']

    return {'RF': RF, 'CO2-fe_C':C, 'CO2-fe_E':E}
