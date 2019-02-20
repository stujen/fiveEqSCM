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
    gas_cycle_parameters.loc['tau1':'tau4'] = np.array([[1000000,394.4,36.54,4.304],[9.,394.4,36.54,4.304],[121.,394.4,36.54,4.304]]).T
    gas_cycle_parameters.loc['r0':'rA'] = np.array([[32.40,0.019,4.165,0.0],\
                  [ 9.05942806e+00, -1.03745809e-07, -1.85711888e-01,  1.45117387e-04],\
                  [ 4.97443512e+01,  5.87120814e-04, -2.02130466e+00,  2.07719812e-02]]).T
    gas_cycle_parameters.loc['PI_conc'] = np.array([278.0,722.0,273.0])
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
              'alpha': pd.DataFrame(alpha.T, index = emissions_in.index, columns = emissions_in.columns)}

    result['RF']['Total'] = result['RF'].sum(axis=1)

    return result

def fit_gas_cycles(emissions_in, \
                   desired_concs, \
                   temp_input, \
                   gas_params = default_gas_params()):


    a = gas_params.loc['a1':'a4'].values.T
    tau = gas_params.loc['tau1':'tau4'].values.T
    r = gas_params.loc['r0':'rA'].values.T
    PI_conc = gas_params.loc['PI_conc'].values
    emis2conc = gas_params.loc['emis2conc'].values
