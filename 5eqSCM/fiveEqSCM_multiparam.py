# 5-equation-SimpleClimateModel code written to include multiple gas inputs, multiple scenario inputs and multiple parameter sets

# Date - 10/02/19
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
import sys
import pandas as pd
# ==================================

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
            return emissions, False
    elif np.array(emissions.shape).size == 1:
        print('Emissions input is 1D (n_timesteps)')
        emissions = emissions[np.newaxis,:]
        emissions = emissions[np.newaxis,...]
    elif np.array(emissions.shape).size == 3 and ((multigas==False or multiscen==False) or (multigas==False and multiscen==False)):
        print('Both flags must be true, since emissions array is 3D')
        return emissions, False
    elif np.array(emissions.shape).size != 3:
        print('Emissions input is not in correct format (needs to be 1D/2D/3D with correct flags)')
        return emissions, False
    
    # check shape of a, r, tau and PI_C parameter arrays are the same size as the number of gases specified
    if (multigas == True) and ((emissions.shape[0] != a.shape[0]) or a.ndim != 2):
        print('multigas = True requires a to have dimensions [num_gases, num_pools], where num_pools = 4')
        print('Currently a has shape: ', a.shape)
        return emissions, False
    elif (multigas == True) and ((emissions.shape[0] != r.shape[0]) or r.ndim != 2):
        print('multigas = True requires r to have dimensions [num_gases, num_pools], where num_pools = 4')
        print('Currently r has shape: ', r.shape)
        return emissions, False
    elif (multigas == True) and ((emissions.shape[0] != tau.shape[0]) or tau.ndim != 2):
        print('multigas = True requires tau to have dimensions [num_gases, num_pools], where num_pools = 4')
        print('Currently tau has shape: ', tau.shape)
        return emissions, False
    elif (multigas == True) and ((emissions.shape[0] != PI_C.shape[0]) or PI_C.ndim != 1):
        print('multigas = True requires PI_C to have dimensions [num_gases]')
        print('Currently PI_C has shape: ', PI_C.shape)
        return emissions, False
    elif (multigas == True) and ((emissions.shape[0] != emis2conc.shape[0]) or emis2conc.ndim != 1):
        print('multigas = True requires emis2conc to have dimensions [num_gases]')
        print('Currently emis2conc has shape: ', emis2conc.shape)
        return emissions, False
    
    return emissions, True

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
    Multiparameterized version of 5eqSCM code
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
    emissions, flag_val = check_format(emissions, multigas, multiscen, a, tau, r, PI_C, emis2conc)
    # if not stop run
    if flag_val == False:
        print('\nRun failed')
        return

    # redefine q if tcr or ecs have changed from default values
    if (tcr != 1.6) * (ecs != 2.75):
        k = 1.0 - (d/70.0)*(1.0 - np.exp(-70.0/d))
        q =  (1.0 / F_2x) * (1.0/(k[0]-k[1])) * np.array([tcr-k[1]*ecs,k[0]*ecs-tcr])
    
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




# ============================================================================================== 

# Define functions to import, run and plot RCP CO2, CH4 and N2O emissions data

def import_RCPs():

    import pandas as pd

    RCP85_E = pd.read_csv('./RCP_data/RCP85_EMISSIONS.csv',skiprows=36,index_col=0)
    RCP85_C = pd.read_csv('./RCP_data/RCP85_MIDYEAR_CONCENTRATIONS.csv',skiprows=37,index_col=0)
    RCP6_E = pd.read_csv('./RCP_data/RCP6_EMISSIONS.csv',skiprows=36,index_col=0)
    RCP6_C = pd.read_csv('./RCP_data/RCP6_MIDYEAR_CONCENTRATIONS.csv',skiprows=37,index_col=0)
    RCP45_E = pd.read_csv('./RCP_data/RCP45_EMISSIONS.csv',skiprows=36,index_col=0)
    RCP45_C = pd.read_csv('./RCP_data/RCP45_MIDYEAR_CONCENTRATIONS.csv',skiprows=37,index_col=0)
    RCP3_E = pd.read_csv('./RCP_data/RCP3PD_EMISSIONS.csv',skiprows=36,index_col=0)
    RCP3_C = pd.read_csv('./RCP_data/RCP3PD_MIDYEAR_CONCENTRATIONS.csv',skiprows=37,index_col=0)

    RCP = {'85':{},'6':{},'45':{},'3':{}}

    RCP['85']['E'] = RCP85_E
    RCP['85']['C'] = RCP85_C
    RCP['6']['E'] = RCP6_E
    RCP['6']['C'] = RCP6_C
    RCP['45']['E'] = RCP45_E
    RCP['45']['C'] = RCP45_C
    RCP['3']['E'] = RCP3_E
    RCP['3']['C'] = RCP3_C

    return RCP

def concplot(C,RCP,rcps):

    import matplotlib
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    from matplotlib import collections
    from matplotlib import pyplot as plt

    plt.rcParams['figure.figsize'] = 16, 9
    plt.rcParams['lines.linewidth'] = 1.5
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.linestyle'] = ':'
    plt.rcParams['grid.color'] = 'black'
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['legend.framealpha'] = 1.0
    plt.rcParams['legend.shadow'] = False
    plt.rcParams['legend.edgecolor'] = 'black'
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['legend.fancybox'] = False

    font = {'weight' : 'normal',
          'size'   : 12}

    plt.rc('font', **font)

    import seaborn as sns

    # -----------------------

    fig,ax=plt.subplots(2,3,figsize=(15,7))
    
    rcps = ['85','6','45','3']
    colors = ['red','brown','green','blue']
    
    for i,rcp in enumerate(rcps):

        ax[0,0].plot(RCP[rcp]['E'].FossilCO2+RCP[rcp]['E'].OtherCO2,color=colors[i],label=rcp+' emissions')
        ax[0,0].set_ylabel('GtC')
        ax[0,0].set_title('CO$_2$')

        ax[1,0].plot(RCP[rcp]['C'].CO2,color=colors[i],label=rcp+' data')
        ax[1,0].plot(RCP[rcp]['C'].CO2.index.values,C[0,i,:],'--',color=colors[i],label=rcp+'-5FAIR')
        ax[1,0].legend(loc='best')
        ax[1,0].set_ylabel('ppm')

        ax[0,1].plot(RCP[rcp]['E'].CH4,color=colors[i])
        ax[0,1].set_ylabel('MtCH$_4$')
        ax[0,1].set_title('CH$_4$')

        ax[1,1].plot(RCP[rcp]['C'].CH4,color=colors[i],label='RCP')
        ax[1,1].plot(RCP[rcp]['C'].CH4.index.values,C[1,i,:],'--',color=colors[i],label='5FAIR')
        ax[1,1].set_ylabel('ppb')

        ax[0,2].plot(RCP[rcp]['E'].N2O,color=colors[i])
        ax[0,2].set_ylabel('MtN$_2$O-N$_2$')
        ax[0,2].set_title('N$_2$O')

        ax[1,2].plot(RCP[rcp]['C'].N2O,color=colors[i],label='RCP')
        ax[1,2].plot(RCP[rcp]['C'].N2O.index.values,C[2,i,:],'--',color=colors[i],label='5FAIR')
        ax[1,2].set_ylabel('ppb')

    plt.tight_layout()

    return fig, ax

def tempplot(T,rcps):

    import matplotlib
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize
    from matplotlib import collections
    from matplotlib import pyplot as plt

    plt.rcParams['figure.figsize'] = 16, 9
    plt.rcParams['lines.linewidth'] = 1.5
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.linestyle'] = ':'
    plt.rcParams['grid.color'] = 'black'
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['legend.framealpha'] = 1.0
    plt.rcParams['legend.shadow'] = False
    plt.rcParams['legend.edgecolor'] = 'black'
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['legend.fancybox'] = False

    font = {'weight' : 'normal',
          'size'   : 12}

    plt.rc('font', **font)

    import seaborn as sns

    # -----------------------

    fig,ax=plt.subplots(figsize=(8,6))
    
    rcps = ['85','6','45','3']
    colors = ['red','brown','green','blue']
    
    for i,rcp in enumerate(rcps):

        ax.plot(np.arange(1765,2501), T[i,:] - np.mean(T[i,1850-1765:1901-1765]),color=colors[i],label=rcp+' U-FaIR temperature response')
        ax.set_ylabel('Temperature anomaly relative to 1850-1900 (K)')
        ax.set_title('T')

    plt.tight_layout()

    return fig, ax

def run_RCPs(rcps = ['85','6','45','3'], plot_out=True):
    '''
    Run RCP emissions scenarios for n gases through U-FAIR simulatneously and plot result
    '''

    import pandas as pd
    from matplotlib import pyplot as plt
    # ----------------------------------

    RCP = import_RCPs()

    emissions = np.zeros((3,4,736))

    for n,rcp_val in enumerate(rcps):
        emissions[0,n,:] = RCP[rcp_val]['E'].FossilCO2.values + RCP[rcp_val]['E'].OtherCO2.values
        emissions[1,n,:] = RCP[rcp_val]['E'].CH4.values
        emissions[2,n,:] = RCP[rcp_val]['E'].N2O.values

    emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12., 16., 28.]) / 28.97)
    a = np.array([[0.2173,0.2240,0.2824,0.2763],[1.,0.,0.,0.],[1.,0.,0.,0.]])
    tau = np.array([[1000000,394.4,36.54,4.304],[9.,394.4,36.54,4.304],[121.,394.4,36.54,4.304]])
    r = np.array([[32.40,0.019,4.165,0.0],\
                 [ 9.05942806e+00, -1.03745809e-07, -1.85711888e-01,  1.45117387e-04],\
                 [ 4.97443512e+01,  5.87120814e-04, -2.02130466e+00,  2.07719812e-02]])
    PI_C = np.array([278.0,722.0,273.0])
    f = np.array([[3.74/np.log(2.),0.,0.],[0,0.,0.036],[0,0,0.12]])

    C,RF,T = multiscen_oxfair(emissions=emissions,emis2conc=emis2conc,a=a,tau=tau,r=r,PI_C=PI_C,f=f,multigas=True,multiscen=True)

    RCP_results = {'C':C, 'RF':RF, 'T':T}

    if plot_out == True:
        fig_C, ax_C = concplot(C,RCP,rcps)
        fig_T, ax_T = tempplot(T,rcps)
        return RCP_results, fig_C, ax_C, fig_T, ax_T
    else:
        return RCP_results











if __name__=='__main__':
    mode = sys.argv[1]

    if mode == 'test':
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib import pyplot as plt
        
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

    if mode == 'simple_multiparam':
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib import pyplot as plt


        fig, ax = plt.subplots()
        ax.set_ylabel('Temperature anomaly relative to 1850-1900 (K)')
        ax.set_title('T')
        ax.set_xlim(1850,2500)
        colors = ['red','brown','green','blue']

        
        param_num = 50 # number of parameter sets 

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

        # import RCP ems scenarios
        rcps = ['85','6','45','3']
        RCP = import_RCPs()
        emissions = np.zeros((3,4,736))
        for n,rcp_val in enumerate(rcps):
            emissions[0,n,:] = RCP[rcp_val]['E'].FossilCO2.values + RCP[rcp_val]['E'].OtherCO2.values
            emissions[1,n,:] = RCP[rcp_val]['E'].CH4.values
            emissions[2,n,:] = RCP[rcp_val]['E'].N2O.values

        for i in return_df.index.values:

            emis2conc = np.array([return_df.emis2conc_CO2.loc[i], return_df.emis2conc_CH4.loc[i], return_df.emis2conc_N2O.loc[i]])
            a = np.array([[return_df.a0_CO2.loc[i],return_df.a1_CO2.loc[i],return_df.a2_CO2.loc[i],return_df.a3_CO2.loc[i]],[return_df.a0_CH4.loc[i],return_df.a1_CH4.loc[i],return_df.a2_CH4.loc[i],return_df.a3_CH4.loc[i]],[return_df.a0_N2O.loc[i],return_df.a1_N2O.loc[i],return_df.a2_N2O.loc[i],return_df.a3_N2O.loc[i]]])
            tau = np.array([[return_df.tau0_CO2.loc[i],return_df.tau1_CO2.loc[i],return_df.tau2_CO2.loc[i],return_df.tau3_CO2.loc[i]],[return_df.tau0_CH4.loc[i],return_df.tau1_CH4.loc[i],return_df.tau2_CH4.loc[i],return_df.tau3_CH4.loc[i]],[return_df.tau0_N2O.loc[i],return_df.tau1_N2O.loc[i],return_df.tau2_N2O.loc[i],return_df.tau3_N2O.loc[i]]])
            r = np.array([[return_df.r0_CO2.loc[i],return_df.rC_CO2.loc[i],return_df.rT_CO2.loc[i],return_df.rA_CO2.loc[i]],\
                          [return_df.r0_CH4.loc[i],return_df.rC_CH4.loc[i],return_df.rT_CH4.loc[i],return_df.rA_CH4.loc[i]],\
                          [return_df.r0_N2O.loc[i],return_df.rC_N2O.loc[i],return_df.rT_N2O.loc[i],return_df.rA_N2O.loc[i]]])
            PI_C = np.array([return_df.PI_C_CO2.loc[i],return_df.PI_C_CH4.loc[i],return_df.PI_C_N2O.loc[i]])
            f = np.array([[return_df.f0_CO2.loc[i],return_df.f1_CO2.loc[i],return_df.f2_CO2.loc[i]],[return_df.f0_CH4.loc[i],return_df.f1_CH4.loc[i],return_df.f2_CH4.loc[i]],[return_df.f0_N2O.loc[i],return_df.f1_N2O.loc[i],return_df.f2_N2O.loc[i]]])
            tcr_value = return_df.TCR.loc[i]
            ecs_value = return_df.ECS.loc[i]
            d = np.array([return_df.d_2.loc[i],return_df.d_1.loc[i]])
            F_2x = return_df.F_2x.loc[i]


            temp_C,temp_RF,temp_T = multiscen_oxfair(emissions=emissions,tcr=tcr_value,ecs=ecs_value,d=d,F_2x=F_2x,emis2conc=emis2conc,a=a,tau=tau,r=r,PI_C=PI_C,f=f,multigas=True,multiscen=True)

            for j,rcp in enumerate(rcps):
                ax.plot(np.arange(1765,2501), temp_T[j,:] - np.mean(temp_T[j,1850-1765:1901-1765]),color=colors[j],label=rcp+' 5eqSCM temperature response', linewidth=1.0)

        fig.savefig('multi_param_test.pdf', dpi=300)
        plt.show()



