import numpy as np
import string
import math
import sys
import pandas as pd

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

from fiveEqSCM_multiparam import *


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

        # fig.savefig('multi_param_test.pdf', dpi=300)
        plt.show()

    if mode == 'compare_to_etminan':
        import scipy as sp
        from scipy.optimize import curve_fit
        def co2_fit_func(co2_conc, a):
            C_0 = 180.0
            return a * np.log(co2_conc / C_0)
        def ch4_fit_func(ch4_conc, a):
            M_0 = 340.0
            return a * (np.sqrt(ch4_conc) - np.sqrt(M_0))
        def n2o_fit_func(n2o_conc, a):
            N_0 = 200.0
            return a * (np.sqrt(n2o_conc) - np.sqrt(N_0))

        fig, ax = plt.subplots(2,2,figsize=(9,7))
        ax[0,0].axhline(y=0, color='black', linewidth=1.5, linestyle='-', alpha=0.5)
        ax[0,1].axhline(y=0, color='black', linewidth=1.5, linestyle='-', alpha=0.5)
        ax[1,0].axhline(y=0, color='black', linewidth=1.5,linestyle='-', alpha=0.5)
        ax[1,1].axhline(y=0, color='black', linewidth=1.5,linestyle='-', alpha=0.5)

        C_0 = 180.0
        N_0 = 200.0
        M_0 = 340.0

        n2o_conc_0 = 323.0
        ch4_conc_0 = 1800.0
        co2_conc_0 = 389.0

        a1 = -2.4e-7
        b1 = 7.2e-4
        c1 = -2.1e-4
        a2 = -8.0e-6
        b2 = 4.2e-6
        c2 = -4.9e-6
        a3 = -1.3e-6
        b3 = -8.2e-6


        # panel a
        co2_conc = np.arange(180,2001)

        etminan_simple_co2_rf = (a1*((co2_conc - C_0)**2) + b1*np.abs(co2_conc - C_0) + c1*(N_0 + n2o_conc_0)/2.0 + 5.36) * np.log(co2_conc / C_0)
        fiveEqSCM_co2_rf = 5.36 * np.log(co2_conc / C_0)

        ax[0,0].plot(co2_conc, fiveEqSCM_co2_rf - fiveEqSCM_co2_rf[co2_conc==co2_conc_0], color='red', label='5eqSCM')
        ax[0,0].plot(co2_conc, etminan_simple_co2_rf - etminan_simple_co2_rf[co2_conc==co2_conc_0], color='black', label='Etminan')
        ax[0,0].set_xlabel('CO$_2$ concentration (ppmv)')
        ax[0,0].set_ylabel('Radiative forcing (Wm$^{-2}$)')
        ax[0,0].text(250,8.85,'(a) CO$_2$')
        # ax[0,0].set_xlim(0,2500)
        # ax[0,0].set_ylim(-6,10)

        p1, p2 = curve_fit(co2_fit_func, co2_conc, etminan_simple_co2_rf)
        ax[0,0].plot(co2_conc, co2_fit_func(co2_conc, p1[0]) - co2_fit_func(co2_conc, p1[0])[co2_conc==co2_conc_0], color='green', label='5eqSCM fit to Etminan', linestyle='-.')
        ax[0,0].legend(loc='lower right', edgecolor='black', framealpha=1.0, fontsize=9)

        # panel b
        ch4_conc = np.arange(340,3501)

        etminan_simple_ch4_rf = (b3*(n2o_conc_0 + N_0)/2.0 + a3*(ch4_conc + M_0)/2.0 + 0.043) * (np.sqrt(ch4_conc) - np.sqrt(M_0))
        fiveEqSCM_ch4_rf = 0.043 * (np.sqrt(ch4_conc) - np.sqrt(M_0))

        ax[0,1].plot(ch4_conc, fiveEqSCM_ch4_rf - fiveEqSCM_ch4_rf[ch4_conc==ch4_conc_0], color = 'red')
        ax[0,1].plot(ch4_conc, etminan_simple_ch4_rf - etminan_simple_ch4_rf[ch4_conc==ch4_conc_0], color = 'black')
        ax[0,1].set_xlabel('CH$_4$ concentration (ppbv)')
        ax[0,1].text(500,0.6,'(b) CH$_4$')
        # ax[0,1].set_xlim(0,3750)
        # ax[0,1].set_ylim(-1.0,0.8)

        q1, q2 = curve_fit(ch4_fit_func, ch4_conc, etminan_simple_ch4_rf)
        ax[0,1].plot(ch4_conc, ch4_fit_func(ch4_conc, q1[0]) - ch4_fit_func(ch4_conc, q1[0])[ch4_conc==ch4_conc_0], color='green', label='5eqSCM fit to Etminan', linestyle='-.')

        # panel c
        n2o_conc = np.arange(200,526)

        etminan_simple_n2o_rf = (a2*(co2_conc_0 + C_0)/2.0 + b2*(n2o_conc + N_0)/2.0 + c2*(ch4_conc_0 + M_0)/2.0 + 0.117) * (np.sqrt(n2o_conc) - np.sqrt(N_0))
        fiveEqSCM_n2o_rf = 0.117 * (np.sqrt(n2o_conc) - np.sqrt(N_0))

        ax[1,0].plot(n2o_conc, fiveEqSCM_n2o_rf - fiveEqSCM_n2o_rf[n2o_conc==n2o_conc_0], color = 'red')
        ax[1,0].plot(n2o_conc, etminan_simple_n2o_rf - etminan_simple_n2o_rf[n2o_conc==n2o_conc_0], color = 'black')
        ax[1,0].set_xlabel('N$_2$O concentration (ppbv)')
        ax[1,0].set_ylabel('Radiative forcing (Wm$^{-2}$)')
        ax[1,0].text(210,0.5,'(c) N$_2$O')
        # ax[1,0].set_xlim(100,600)
        # ax[1,0].set_ylim(-0.45,0.65)

        r1, r2 = curve_fit(n2o_fit_func, n2o_conc, etminan_simple_n2o_rf)
        ax[1,0].plot(n2o_conc, n2o_fit_func(n2o_conc, r1[0]) - n2o_fit_func(n2o_conc, r1[0])[n2o_conc==n2o_conc_0], color='green', label='5eqSCM fit to Etminan', linestyle='-.')

        # panel d
        n2o_conc_0 = 525.0
        ch4_conc_0 = 3500.0
        co2_conc_0 = 2000.0

        n2o_conc_1 = 323.0
        ch4_conc_1 = 1800.0
        co2_conc_1 = 389.0


        n2o_overlap_of_co2_rf = c1 * (n2o_conc_0 - N_0) * np.log(co2_conc/C_0)
        n2o_overlap_of_ch4_rf = b3 * (n2o_conc_0 - N_0) * (np.sqrt(ch4_conc) - np.sqrt(M_0))
        ch4_overlap_of_n2o_rf = c2 * (ch4_conc_0 - M_0) * (np.sqrt(n2o_conc) - np.sqrt(N_0))
        co2_overlap_of_n2o_rf = a2 * (co2_conc_0 - C_0) * (np.sqrt(n2o_conc) - np.sqrt(N_0))

        ax[1,1].plot(co2_conc, n2o_overlap_of_co2_rf - n2o_overlap_of_co2_rf[co2_conc==co2_conc_1], color = 'black', linestyle='--', label='N$_2$O overlap of CO$_2$')
        ax[1,1].plot(ch4_conc, n2o_overlap_of_ch4_rf - n2o_overlap_of_ch4_rf[ch4_conc==ch4_conc_1], color = 'red', linestyle='--', label='N$_2$O overlap of CH$_4$')
        ax[1,1].plot(n2o_conc, ch4_overlap_of_n2o_rf - ch4_overlap_of_n2o_rf[n2o_conc==n2o_conc_1], color = 'green', linestyle=':', label='CH$_4$ overlap of N$_2$O')
        ax[1,1].plot(n2o_conc[::50], co2_overlap_of_n2o_rf[::50] - co2_overlap_of_n2o_rf[n2o_conc==n2o_conc_1], color = 'green', marker='+', linestyle='', label='CO$_2$ overlap of N$_2$O')
        ax[1,1].set_ylim(-0.35,0.1)
        ax[1,1].set_xlim(0,3750)
        ax[1,1].text(1450,0.05,'(d) overlap strength')
        ax[1,1].set_xlabel('Overlapped gas concentration (ppmv/ppbv)')
        ax[1,1].legend(loc='lower left', edgecolor='black', framealpha=1)

        # fig.savefig('etminan_reproduction_5eqSCM.pdf',dpi=300)
        plt.show()

        # with plt.xkcd():
        #     # Based on "Stove Ownership" from XKCD by Randall Monroe
        #     # http://xkcd.com/418/

        #     fig = plt.figure()
        #     ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
        #     ax.spines['right'].set_color('none')
        #     ax.spines['top'].set_color('none')
        #     plt.xticks([])
        #     plt.yticks([])
        #     ax.set_ylim([-30, 10])

        #     data = np.ones(100)
        #     data[70:] -= np.arange(30)

        #     plt.annotate(
        #         'THE DAY I REALIZED\nI COULD COOK BACON\nWHENEVER I WANTED',
        #         xy=(70, 1), arrowprops=dict(arrowstyle='->'), xytext=(15, -10))

        #     plt.plot(data)

        #     plt.xlabel('time')
        #     plt.ylabel('my overall health')
        #     fig.text(
        #         0.5, 0.05,
        #         '"Stove Ownership" from xkcd by Randall Monroe',
        #         ha='center')

        #     plt.show()

    if mode == 'check_HC134a':

        import scipy as sp

        # import RCP ems scenarios
        rcps = ['85','6','45','3']
        RCP = import_RCPs()
        emissions = np.zeros((4,736))

        for n,rcp_val in enumerate(rcps):
            emissions[n,:] = RCP[rcp_val]['E'].HFC134a.values

        
        ## Testing minimize based fitter (seems to work):
        def rcpdiff(x, rcps):

            diff = 0

            emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.,16.,28.]) / 28.97),
            a = np.array([[0.2173,0.2240,0.2824,0.2763],[1.,0.,0.,0.],[1.,0.,0.,0.]]),
            tau = np.array([[1000000,394.4,36.54,4.304],[9.,394.4,36.54,4.304],[121.,394.4,36.54,4.304]]),
            PI_C = np.array([278.0,722.0,273.0]),
            iirf100_max = 97.0,
            f = np.array([[3.74/np.log(2.),0.,0.],[0,0.,0.036],[0,0,0.12]])
            
            for i, rcp_test in enumerate(rcps):
                # should we be fitting in isolation????
                C = multiscen_oxfair(emissions=emissions[i,:], r=np.array([[x[0],x[1],x[2],x[3]]]), multiscen=False, multigas=False)[0]
                
                print(C)
                # diff += np.sum((C - RCP[rcp_test]['C'].CH4.values)**2)
                
            return diff

        rcpdiff((10,0.0,0.0,0.0),['3'])
        # sp.optimize.minimize(rcpdiff,(10,0.0,0.0,0.0))

