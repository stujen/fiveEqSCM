# Universal-FaIR code written so it can take multiple gas inputs and multiple scenario inputs

# Date - 30/12/18
# Written by Stuart Jenkins (stuart.jenkins@wadham.ox.ac.uk) and Nicholas Leach (nicholas.leach@stx.ox.ac.uk) 

# ===========================================================================================================

# TO DO...
# define function to unpack pandas input into formats required?
# define n_gases, n_scens, n_timesteps and use throughout?
# check default inputs are correct?



# import required external packages
import numpy as np
import string
import math

# ==================================

def check_format(emissions, multigas, multiscen, a, tau, r, PI_C, emis2conc):
    
    # check if emissions is 3D, if not make it 3D, complying with shape requirements
    if np.array(emissions.shape).size == 2:
        if multigas == True and multiscen == False:
            print('Emissions input is 2D (n_gases, n_timesteps)')
            emissions = emissions[:,np.newaxis,:]
        elif multigas == False and multiscen == True:
            print('Emissions input is 2D (n_scenarios, n_timesteps)')
            emissions = emissions[np.newaxis,...]
        elif (multigas == False and multiscen == False) or (multigas == True and multiscen == True):
            print('One flag (multigas/multiscen) must be switched True, other False')
            return False
    elif np.array(emissions.shape).size == 1:
        print('Emissions input is 1D (n_timesteps)')
        emissions = emissions[np.newaxis,:]
        emissions = emissions[np.newaxis,...]
    elif np.array(emissions.shape).size == 3 and ((multigas==False or multiscen==False) or (multigas==False and multiscen==False)):
        print('Both flags must be true, since emissions array is 3D')
        return False
    elif np.array(emissions.shape).size != 3:
        print('Emissions input is not in correct format (needs to be 1D/2D/3D with correct flags)')
        return False
    
    # check shape of a, r, tau and PI_C parameter arrays are the same size as the number of gases specified
    if (multigas == True) and ((emissions.shape[0] != a.shape[0]) or a.ndim != 2):
        print('multigas = True requires a to have dimensions [num_gases, num_pools], where num_pools = 4')
        print('Currently a has shape: ', a.shape)
        return False
    elif (multigas == True) and ((emissions.shape[0] != r.shape[0]) or r.ndim != 2):
        print('multigas = True requires r to have dimensions [num_gases, num_pools], where num_pools = 4')
        print('Currently r has shape: ', r.shape)
        return False
    elif (multigas == True) and ((emissions.shape[0] != tau.shape[0]) or tau.ndim != 2):
        print('multigas = True requires tau to have dimensions [num_gases, num_pools], where num_pools = 4')
        print('Currently tau has shape: ', tau.shape)
        return False
    elif (multigas == True) and ((emissions.shape[0] != PI_C.shape[0]) or PI_C.ndim != 1):
        print('multigas = True requires PI_C to have dimensions [num_gases]')
        print('Currently PI_C has shape: ', PI_C.shape)
        return False
    elif (multigas == True) and ((emissions.shape[0] != emis2conc.shape[0]) or emis2conc.ndim != 1):
        print('multigas = True requires emis2conc to have dimensions [num_gases]')
        print('Currently emis2conc has shape: ', emis2conc.shape)
        return False
    
    return emissions

def g_1(a,tau,h):
    
    g1 = np.sum( a*tau*(1. - (1.+h/tau)*np.exp(-h/tau)), axis=1)
    
    return g1

def g_0(a,tau,h):
    
    g0 = (np.sinh( np.sum( a*tau*(1. - np.exp(-h/tau)), axis=1) / g_1(a,tau,h)) )**(-1.)
    
    return g0

def alpha_val(G,G_A,T,tau,a,r,pre_ind_C,h=100,iirf100_max=97.0):
    
    iirf100_val = r[...,np.newaxis,0] + r[...,np.newaxis,1]*(G-G_A) + r[...,np.newaxis,2]*T + r[...,np.newaxis,3]*G_A
    
    # if iirf100 value larger than max value, set to max value 
    iirf100_val = (iirf100_val>iirf100_max)*iirf100_max + iirf100_val*(iirf100_val<iirf100_max)
    
    alpha_val = g_0(a,tau,h)[...,np.newaxis] * np.sinh(iirf100_val / g_1(a,tau,h)[...,np.newaxis])
    
    return alpha_val

def step_conc(R,alpha,E,a,tau,pre_ind_C,emis2conc):
    
    R = E[...,np.newaxis] * emis2conc[:,np.newaxis,np.newaxis]*a[:,np.newaxis,:]*alpha[...,np.newaxis]*tau[:,np.newaxis,:] * (1. - np.exp(-1./(alpha[...,np.newaxis]*tau[:,np.newaxis,:]))) + R*np.exp(-1./(alpha[...,np.newaxis]*tau[:,np.newaxis,:]))
    
    C = pre_ind_C[:,np.newaxis] + np.sum(R,axis=2)
    
    G_A = (C - pre_ind_C[:,np.newaxis]) / emis2conc[:,np.newaxis]
    
    return C,R,G_A

def step_forc(C,pre_ind_C,F_ext,f):
    
    F = np.sum(f[...,np.newaxis,0]*np.log(C/pre_ind_C[:,np.newaxis]) + f[...,np.newaxis,1]*(C - pre_ind_C[:,np.newaxis]) + f[...,np.newaxis,2] * (np.sqrt(C) - np.sqrt(pre_ind_C[:,np.newaxis])), axis=0) + F_ext
    
    return F

def step_temp(S,F,q,d):
    
    S = q[np.newaxis,:]*F[:,np.newaxis]*(1-np.exp(-1/d[np.newaxis,:])) + S*np.exp(-1/d[np.newaxis,:])
    
    T = np.sum(S, axis=1)
    
    return S,T

def multiscen_oxfair(emissions,
                     emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.,16.,28.]) / 28.97),
                     a = np.array([[0.2173,0.2240,0.2824,0.2763],[1.,0.,0.,0.],[1.,0.,0.,0.]]),
                     tau = np.array([[1000000,394.4,36.54,4.304],[9.,394.4,36.54,4.304],[121.,394.4,36.54,4.304]]),
                     r = np.array([[32.40,0.019,4.165,0.0],\
                                   [ 9.05942806e+00, -1.03745809e-07, -1.85711888e-01,  1.45117387e-04],\
                                   [ 4.97443512e+01,  5.87120814e-04, -2.02130466e+00,  2.07719812e-02]]),
                     PI_C = np.array([278.0,722.0,273.0]),
                     iirf100_max = 97.0,
                     f = np.array([[3.74/np.log(2.),0.,0.],[0,0.,0.036],[0,0,0.12]]),
                     tcr = 1.6,
                     ecs = 2.75,
                     d = np.array([239.0,4.1]),
                     q = np.array([0.33,0.41]),
                     F_2x = 3.74,
                     multigas=False, 
                     multiscen=False):
    '''
    Multiparameterized version of (OxFaIR/U-FaIR) code.
    ===================================================

    Inputs
    ------

    emissions : emissions take shape 1D timeseries, 2D [n_gases, n_timesteps] or 2D [n_scenarios, n_timesteps], 3D [n_gases, n_scenarios, n_timesteps].
                If 2D emissions array is inputted user must switch relevent flags to True so U-FaIR knows how to treat array (multigas/multiscen).
                For 3D emissions array user must switch multigas and multiscen to True.
                In all cases correct parameter array shapes must be supplied with the given emissions shape.

    Additional parameters
    ---------------------

    emis2conc : conversion factor between emissions and concentrations units. Size = [n_gases]

    a :         pool splits for concentrations. Maximum is 4 pool model for gases. 
                CH4 and N2O have one pool model [1.,0.,0.,0.], CO2 has 4 pool model [0.2173,0.2240,0.2824,0.2763].
                Size = [n_gases, 4]

    tau :       Time constants for the decay rates of the 4 pools with splits given by a. 
                Size = [n_gases, 4].

    r :         iIRF100 estimation parameter values (r0, rC, rT, rA). e.g. standard values for CO2 are [32.4, 0.019, 4.165, 0.0]
                Size = [n_gases, 4].

    PI_C :      Pre-Industrial Concentrations (PI_C) for the n_gases in model. 
                Size = [n_gases].

    iirf100_max : Max allowed value of integrated impulse response function (h=100yrs) before saturation. Default = 97.0

    f :         Contributions from log, linear and sqrt forcingfor each gas in n_gases. 
                Size = [n_gases, 3].

    tcr :       TCR value for run. (Default = 1.6).

    ecs :       ECS value for run. (Default = 2.75).

    d :         2 box temperature model time constants. (Default = [239.,4.1]).

    q :         Alternate representation of TCR and ECS values. Fractional split between temperature boxes. (Default = []0.33, 0.41).

    F_2x :      Forcing change for 2x CO2 concentration change. (Default = 3.74).

    multigas :  Flag for multigas run. True if run includes multiple gases emissions, False otherwise. (Default = False).

    multiscen :  Flag for multiscen run. True if run includes multiple emissions scnearios, False otherwise. (Default = False).

    '''
    
    # ------------------
    
    # check format of inputs is correct
    emissions = check_format(emissions, multigas, multiscen, a, tau, r, PI_C, emis2conc)
    # if not stop run
    if emissions == False:
        print('\nRun failed')
        return
    
    # define empty arrays to populate
    # -------------------
    G = np.cumsum(emissions,axis=2) # cumulative emissions array
    C = np.zeros(emissions.shape) # concentrations array
####### currently calculate a total RF for all gases in each scenario together, disgree with the total RF step
    RF = np.zeros(emissions[0].shape) # RF array
#######
    T = np.zeros(emissions[0].shape) # total temperature response array
    alpha = np.zeros(emissions.shape) # alpha array
    
    # calculate zeroth step values
    # -------------------
    # calculate alpha values for all gases and all scenarios for zeroth timestep
    alpha[...,0] = alpha_val(G=np.zeros((emissions.shape[0],emissions.shape[1])),
                            G_A=np.zeros((emissions.shape[0],emissions.shape[1])),
                            T=np.zeros((emissions.shape[0],emissions.shape[1])),
                            tau=tau,a=a,r=r,h=100.,pre_ind_C=PI_C,iirf100_max = 97.0)
    # calculate concentrations (C), concentration split betwenn pools (R), 
    #    and cumulative emissions remaining in atmosphere (G_A) for zeroth timestep
    C[...,0],R,G_A = step_conc(R=np.zeros((emissions.shape[0],emissions.shape[1],4)),
                            alpha=alpha[...,0],E=emissions[...,0],a=a,tau=tau,
                            pre_ind_C=PI_C,emis2conc=emis2conc)
    # calculate RF array for zeroth timestep
    RF[...,0] = step_forc(C=C[...,0],pre_ind_C=PI_C,F_ext=0.,f=f)
    # calculate temperature response (T) and split between boxes (S) for zeroth timestep
    S,T[...,0] = step_temp(S=np.zeros((emissions.shape[1],2)),F=RF[...,0],q=q,d=d)
    
    # run over all other timesteps...
    for t in np.arange(1,emissions.shape[2]):
        
        # calculate alpha values for all gases and all scenarios for t'th timestep
        alpha[...,t] = alpha_val(G=G[...,t-1],G_A=G_A,T=T[...,t-1],tau=tau,a=a,r=r,h=100.,pre_ind_C=PI_C,iirf100_max = 97.0)
        # calculate concentrations (C), concentration split betwenn pools (R), 
        #    and cumulative emissions remaining in atmosphere (G_A) for t'th timestep
        C[...,t],R,G_A = step_conc(R = R,alpha=alpha[...,t],E=emissions[...,t],a=a,tau=tau,pre_ind_C=PI_C,emis2conc=emis2conc)
        # calculate RF array for t'th timestep
        RF[...,t] = step_forc(C=C[...,t],pre_ind_C=PI_C,F_ext=0.,f=f)
        # calculate temperature response (T) and split between boxes (S) for t'th timestep
        S,T[...,t] = step_temp(S=S,F=RF[...,t],q=q,d=d)
        
    # return C, RF and T arrays
    return C,RF,T

