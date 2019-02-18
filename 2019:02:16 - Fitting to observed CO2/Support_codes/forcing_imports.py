## File contains functions to import forcing datasets, sample TCR/ECS parameter distributions, calculate CC parameters in fair runs

# Written by Stuart Jenkins (stuart.jenkins@wadham.ox.ac.uk)
# -----------------------------------------------------------------------------------

#Imports
# -------------------------------------------------
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.special import erf, erfinv
from scipy.optimize import root
from scipy.stats import norm, beta
from statsmodels.api import OLS
import statsmodels.tools.tools
import datetime
import os
import h5py
import copy
from pandas import DataFrame
# -------------------------------------------

# Import the FaIR climate-carbon-cycle model from python file
from fair_scm import *
# -------------------------------------------------

# Definitions
# -------------------------------------------------
# Set shape of TCR CDF
shape='gaussian'
# set key parameters
a=3.74/np.log(2.0)
c0=278.
# set baseline period for FaIR model runs
base_low=1861.
base_high=1880.
# Set FaIR model temperature boxes parameters
d=np.array([409.5,8.4])
years = np.arange(1765.0,2501.0)
# -------------------------------------------------





# -----------------------------------------------
# return the timeseries x re-baselined to be zero between the years base_low and base_high.
# -----------------------------------------------
def express_as_anom(x,years,base_low=1861.,base_high=1880.):
    # Inputs are a timeseries (x) and corresponding years, and a baseline period to rebaseline the data to
    # Returns the timeseries (x) rebaselined to the specified reference period
    return x - np.mean(x[np.logical_and(years>=base_low,years<=base_high)])

# -------------------------------------------------------------------- #


# -----------------------------------------------
# take the rough bounds of a gaussians central 2/3rds range, and return the corresponding gaussians mean and S.D.
# -----------------------------------------------
def analy_musig(bound1,bound2):
    # Inputs the values of a gaussian distribution at the boundaries of the (central 2/3rds) likely range
    # Calculates the mean and standard deviation of the gaussian which fits these bounds

    C = (np.sqrt(2.0)*erfinv((0.5-1.0/3.0)*2 - 1.0 ))
    K = erfinv((0.5+1.0/3.0)*2 - 1)*np.sqrt(2.0)

    mu_ans = ( K* np.log(bound1) - np.log(bound2))/(K/C -1) 
    sig_ans = (np.log(bound1) - mu_ans) / C

    mu_ans = ( np.log(bound1)*erfinv(2*(0.5+1.0/3.0)-1.0) - np.log(bound2)*erfinv(2*(0.5-1.0/3.0)-1.0)  ) / (erfinv(2*(0.5+1.0/3.0)-1.0) - erfinv(2*(0.5-1.0/3.0)-1.0))
    sig_ans = ( np.log(bound2) - mu_ans) / (np.sqrt(2.0)*erfinv(2*(0.5+1.0/3.0)-1.0)  )

    return mu_ans, sig_ans

# -------------------------------------------------------------------- #


# -----------------------------------------------
# load the historical contrails radiative forcing dataset from AR5 data file.
# -----------------------------------------------
def load_contrails_rf(years,datadir='./Data/'):
    
    rf_cont_file = datadir+'ipcc_tod_rftimeseries.txt'
    hist_rf_dat = np.genfromtxt(rf_cont_file,skip_header=4)
    #Contrail forcing assumed to be constant in future as no scenario data 
    rf_cont = 0.05*np.ones_like(years)
    rf_cont[years<=2011] = hist_rf_dat[15:,10]

    return rf_cont

# -------------------------------------------------------------------- #


# -----------------------------------------------
# load the components of the RCP8.5 radiative forcing dataset, 
#     return datasets in arrays of natural forcing contributions and anthropogenic forcing contributions.
# -----------------------------------------------
def load_rf_comps(rf_file,ar5_scale=False):
    """Loads the components of RF from a spreadsheet (.csv) of the different components
    """
    
    # outline the columns layout in the rf data file
    dt = np.dtype({'names':["YEARS","TOTAL_INCLVOLCANIC_RF","VOLCANIC_ANNUAL_RF","SOLAR_RF","TOTAL_ANTHRO_RF","GHG_RF",
    "KYOTOGHG_RF","CO2CH4N2O_RF","CO2_RF","CH4_RF","N2O_RF","FGASSUM_RF","MHALOSUM_RF","CF4","C2F6","C6F14","HFC23","HFC32",
    "HFC43_10","HFC125","HFC134a","HFC143a","HFC227ea","HFC245fa","SF6","CFC_11","CFC_12","CFC_113","CFC_114","CFC_115",
    "CARB_TET","MCF","HCFC_22","HCFC_141B","HCFC_142B","HALON1211","HALON1202","HALON1301","HALON2402","CH3BR","CH3CL",
    "TOTAER_DIR_RF","OCI_RF","BCI_RF","SOXI_RF","NOXI_RF","BIOMASSAER_RF","MINERALDUST_RF","CLOUD_TOT_RF","STRATOZ_RF",
    "TROPOZ_RF","CH4OXSTRATH2O_RF","LANDUSE_RF","BCSNOW_RF"],'formats':54*["f8"]})
    #load foricng data from file
    forc_data = np.genfromtxt(rf_file,skip_header=59,delimiter=',',dtype=dt)
    
    # split into individual components 
    years = forc_data["YEARS"]
    rf_co2 = forc_data["CO2_RF"]
    rf_wmghg = forc_data['GHG_RF']
    rf_owmghg = rf_wmghg - rf_co2
    rf_landuse = forc_data["LANDUSE_RF"]
    rf_bcos = forc_data["BCSNOW_RF"]
    rf_oz = (forc_data["TROPOZ_RF"]+forc_data["STRATOZ_RF"])
    rf_strath2o = forc_data["CH4OXSTRATH2O_RF"]
    rf_aero = (forc_data["TOTAER_DIR_RF"] + forc_data["CLOUD_TOT_RF"])

    # output array of anthro-rf contributions
    rf_out = np.array([rf_co2,rf_owmghg,rf_oz,rf_strath2o,rf_landuse,rf_bcos,rf_aero])
    
    rf_solar = forc_data["SOLAR_RF"]
    #Smooth future solar cycle forcing with a ten year gaussian filter
    rf_solar[years>2015] = gaussian_filter1d(forc_data["SOLAR_RF"][years>2015],10)
    
    # output array of natural-rf contributions
    rf_nat = np.array([rf_solar,forc_data["VOLCANIC_ANNUAL_RF"]])
    
    # if AR5 scale is set True scale the anthro forcing contributions so they match AR5 quoted values in 2011
    if ar5_scale == True:
        #Scale the anthropgogenic forcing to get the AR5 values in 2011
        # Entry 1 in list is the total ERF for WMGHG
        rf_2011_mid = np.array([1.82,2.83,0.35,0.07,-0.15,0.04,-0.9])
        #Set the first component to be for other WMGHG (not including CO2) only
        rf_2011_mid[1] = rf_2011_mid[1] - rf_2011_mid[0]
        #Find the set of scaling factors needed to get AR5 forcings in 2011
        req_sf = rf_2011_mid[:,np.newaxis] / rf_out[:,years==2011]
        rf_out = req_sf * rf_out
        
        #Scale the solar forcing to AR5 values, leave volcanic as RCP timeseries
        rf_nat[0] = 0.05 / rf_nat[0,years==2011] * rf_nat[0]
        return rf_out, rf_nat, req_sf
    else:
        # else don't scale them to AR5 values and just output result
        return rf_out, rf_nat
    
    #Add on the AR5 forcings before 2011 (from 'ipcc_tod_rftimeseries.txt' dataset)
    dt = np.dtype({'names':["YEARS","CO2_RF","OWMGHG_RF","TROPOZ_RF","STRATOZ_RF","TOTAER_DIR_RF","TOTAER_RF","LANDUSE_RF",
    "CH4OXSTRATH2O_RF","BCSNOW_RF","CONTRAIL_RF","SOLAR_RF","VOLCANIC_ANNUAL_RF"],'formats':13*["f8"]})
    forc_data = np.genfromtxt('./Data/ipcc_tod_rftimeseries.txt',skip_header=4,dtype=dt)

    ar5_years = forc_data["YEARS"][15:]
    rf_co2 = forc_data["CO2_RF"][15:]
    rf_owmghg = forc_data['OWMGHG_RF'][15:]
    rf_landuse = forc_data["LANDUSE_RF"][15:]
    rf_bcos = forc_data["BCSNOW_RF"][15:]
    rf_oz = (forc_data["TROPOZ_RF"][15:]+forc_data["STRATOZ_RF"][15:])
    rf_strath2o = forc_data["CH4OXSTRATH2O_RF"][15:]
    rf_aero = forc_data["TOTAER_RF"][15:]

    rf_out_ar5 = np.array([rf_co2,rf_owmghg,rf_oz,rf_strath2o,rf_landuse,rf_bcos,rf_aero])

    rf_out[:,years<=2011] = rf_out_ar5

    if ar5_scale == True:
        return rf_out, rf_nat, req_sf
    else:
        return rf_out, rf_nat
    
# -------------------------------------------------------------------- #


# -----------------------------------------------
# combine the non-CO2 gaussian forcing contributions into one timeseries 'rf_gauss'
# -----------------------------------------------
def comb_rf_comps(rf_scens,rf_comps,name):
    #Sum the "guassian" components of the forcing together
    rf_scens[name]={}
    rf_scens[name]['rf_co2'] = rf_comps[0,:]
    rf_scens[name]['rf_gauss'] = np.sum(rf_comps[1:-2,:],axis=0)
    rf_scens[name]['rf_aero'] = rf_comps[-1,:]
    rf_scens[name]['rf_bc'] = rf_comps[-2,:]
    
    return rf_scens

# -------------------------------------------------------------------- #


# -----------------------------------------------
# create the distributions for the TCR and ECS values, and then sample them for the required CDF values (mean and central 1/3rd and 2/3rds)
# -----------------------------------------------
def create_climresp_dist(tcr_lr=[1.0,2.5],ecs_lr=[1.5,4.5]):
    """ Create the sampling of the climate response
        parameter distributions """

    #Set the points in the cumulative density function to sample
    cdfs = np.array([0.5-(1.0/3.0),0.5,0.5+(1.0/3.0),1.0/3.0,2.0/3.0])

    if shape == 'lognorm':

        #Fit the log-normal distributions
        mu_tcr, sig_tcr = analy_musig(tcr_lr[0],tcr_lr[1])
        mu_ecs, sig_ecs = analy_musig(ecs_lr[0],ecs_lr[1])

        #Compute the TCR and ECS intervals that are consistent 
        #with likely above and likely below          
        TCR = np.exp( sig_tcr * np.sqrt(2.0) * erfinv(2*(cdfs) -1.0) + mu_tcr)
        ECS = np.exp( sig_ecs * np.sqrt(2.0) * erfinv(2*(cdfs) -1.0) + mu_ecs)


    #Gaussian TCR  distributions
    if shape == 'gaussian':
        #Use symmettry to find the mean of the TCR distribution
        mu_tcr = 0.5 * (tcr_lr[0] + tcr_lr[1])
        TCR = [1.0,1.75,2.5]
        #Fit lognormal ECS distribution 
        mu_ecs, sig_ecs = analy_musig(1.5,4.5)
        #Find the TCR Standard Deviation 
        sig_tcr = (1.0 - 1.75)/(erfinv((0.5-(1.0/3.0))*2 - 1)*np.sqrt(2.0))
        #Compute the TCR and ECS intervals that are consistent 
        #with likely above and likely below
        TCR = sig_tcr * np.sqrt(2.0) * erfinv(2*(cdfs) -1.0) + mu_tcr
        ECS = np.exp( sig_ecs * np.sqrt(2.0) * erfinv(2*(cdfs) -1.0) + mu_ecs)

    return cdfs, TCR, ECS

# -------------------------------------------------------------------- #


# -----------------------------------------------
# calculate the parameter sets required for the carbon cycle, run with each scenario 
#     with the different TCR/ECS values. Parameters sets are found so that the correct 
#     CO2 emissions in 2017 are found.
# -----------------------------------------------
def calc_cc_uncert(TCR, ECS, cdfs, rf_scens, sf_aero, sf_gauss, rf_cont, rf_nat, name='RCP8.5',start_r0=35.0,start_rT=4.5,start_rC=0.02):
    """ Calculates the adaptive blended temperature pathways to get to a 2100 goal
    and their compatible CO2 emissions """

    #Load the data from sampling carbon cylce uncertainty 
    #Create gaussian that spans CMIP5 rC:r0 ratio with central value of 1...
    mu_cc = 1.0
    sig_cc = (0.0 - 1.0)/(erfinv((0.5-(1.0/3.0))*2 - 1)*np.sqrt(2.0))
    #Sample this distribution for all the cdfs considered... 
    cc_ratfeed = sig_cc * np.sqrt(2.0) * erfinv(2*np.array(cdfs) -1.0) + mu_cc

    # create arrays to store output parameter values
    s_temp = np.zeros_like(TCR)
    s_cu = np.zeros_like(TCR)
    iirf_100_preind = np.zeros_like(TCR)

    for i in range(0,len(TCR)):
            # for each TCR and ECS value, scale the forcing appropriately to get total non-CO2 forcing
            nco2_forc = sf_gauss[i]*rf_scens[name]['rf_gauss'] +\
                        rf_scens[name]['rf_bc'] + sf_aero[i]*rf_scens[name]['rf_aero'] +\
                        rf_cont + rf_nat
            #Find carbon cycle feedback scaling to get 2017 emissions correct...
            funt = lambda s: fair_scm_emsback(rf_scens[name]['rf_co2'],other_rf= nco2_forc,TCR=float(TCR[i]),ECS=float(ECS[i]),a=a,
                                              rT=s*cc_ratfeed[i]* start_rT,r0=s*start_r0,rC = s*cc_ratfeed[i]*start_rC,
                                              C_0=c0,d1=d[1],d2=d[0])[1][years==2017] - 11.3
            s_o = root(funt,1.0).x

            #Save the carbon cycle parameters 
            s_temp[i] = s_o * cc_ratfeed[i] *  start_rT
            s_cu[i] = s_o * cc_ratfeed[i] *  start_rC
            iirf_100_preind[i] = s_o *  start_r0

    return s_temp, s_cu, iirf_100_preind

# -------------------------------------------------------------------- #


# -----------------------------------------------
# calculate the TCRE value for a given TCR/ECS value
# -----------------------------------------------
def calc_tcre(TCR, ECS, iirf_100_preind, s_cu, s_temp):
    """Calculate the TCRE from the parameter sets used"""
    #Set up the 1%/yr concentration ramp
    t = np.arange(0,600.0)

    conc_opc = c0 * np.exp( np.log(1.01) * (t))
    rf_opc = a * np.log(conc_opc/c0)

    TCRE = np.zeros_like(TCR)

    for i in range(0,len(TCR)):
        temps_opc, emms_opc = fair_scm_emsback(rf_opc,other_rf= np.zeros_like(rf_opc),TCR=float(TCR[i]),
                ECS=float(ECS[i]),a=a,rT=float(s_temp[i]),
                r0=float(iirf_100_preind[i]),rC = float(s_cu[i]),
                C_0=c0,d1=d[1],d2=d[0])
        cems_opc = np.cumsum(emms_opc)

        TCRE[i] = temps_opc[np.argmin(np.abs(cems_opc-1000.0))]

    return TCRE
# -------------------------------------------------------------------- #


# -----------------------------------------------
# print the important parameters for the given TCR/ECS values used
# -----------------------------------------------
def print_parameters(cdfs, TCR, ECS, TCRE, rf_aero_2011, rf_tot_2011, sf_gauss, sf_aero, iirf_100_preind, s_cu, s_temp):
    
    f_out = open('./Output/parameters_out.txt','w')
    f_out.write('Parameters for the finalised SPM figure 1 code given to IPCC'+'\n\n')
    f_out.write('cdf'+'\t'+'TCR'+'\t'+'ECS'+'\t'+'TCRE'+'\t'+'rf_aero_2011'+'\t'+'rf_tot_2011'+'\t'+'sf_gauss'+'\t'+'sf_aero'+'\t'+'r0'+'\t'+'rC'+'\t'+'rT'+'\n')

    for i in range(0,len(TCR)):
        f_out.write(str(cdfs[i])+'\t'+str(TCR[i])+'\t'+str(ECS[i])+'\t'+str(TCRE[i])+'\t'+str(rf_aero_2011[i])+'\t'+str(rf_tot_2011[i])+'\t'+str(sf_gauss[i])+'\t'+str(sf_aero[i])+'\t'+str(iirf_100_preind[i])+'\t'+str(s_cu[i])+'\t'+str(s_temp[i])+'\n')
    f_out.close()
    
    print 'CDFs:\t\t ', np.round(cdfs,2)
    print 'TCRs:\t\t ', np.round(TCR,2)
    print 'ECSs:\t\t ', np.round(ECS,2)
    print 'TCREs:\t\t ', np.round(TCRE,2)
    print 'rf_aero_2011:\t ', np.round(rf_aero_2011,2)
    print 'rf_tot_2011:\t ', np.round(rf_tot_2011,2)
    print 'sf_gauss:\t ', np.round(sf_gauss,2)
    print 'sf_aero:\t ', np.round(sf_aero,2)
    print 'r0:\t\t ', np.round(iirf_100_preind,2)
    print 'rC:\t\t ', np.round(s_cu,3)
    print 'rT:\t\t ', np.round(s_temp,3)
        
    return












# -----------------------------------------------
# Scaling of gaussian non-CO2 forcing component
# Scale gaussian non-CO2 forcing to span the 5-95% certainty range quoted in AR5 (assuming +/-20% uncertainty on WMGHG forcings)
# Use the TCR CDF shape to define the gaussian scaling factor CDF shape to sample
# Methodology from Millar et al. (Nature Geoscience, 2017)
# -----------------------------------------------
def create_gauss_scalingfactors(cdfs):
    """ Create the scaling factor distributions needed to sample
        uncertainty in the Gaussian non-CO2 radiative forcings"""
    #Scale based on combined gaussian components 
    rf_2011_mid = np.array([1.82,2.83,0.35,0.07,-0.15,0.04,-0.9])
    rf_2011_up = [2.18,3.4,0.559,0.121,-0.047,0.09,-0.1]
    rf_2011_low = [1.46,2.26,0.141,0.019,-0.253,0.019,-1.9]
    #Estimate sd using 5-95 intervals
    erf_sigs = (np.array(rf_2011_up) - np.array(rf_2011_low)) / (2*1.654)
    sig_wmghg = np.copy(erf_sigs[1])
    #Find the non-CO2 GHG forcing uncertainty 
    sig_owmghg = np.sqrt(erf_sigs[1]**2 - erf_sigs[0]**2)
    erf_sigs[1] = sig_owmghg
    sig_tot = np.sqrt(np.sum(erf_sigs[1:-2]**2))
    rf_2011_mid_a = np.copy(rf_2011_mid)
    rf_2011_mid_a[1] = rf_2011_mid[1] - rf_2011_mid[0]
    #Calculate the scaling factors
    #Derive the scaling factors to span 5-95% AR5 guassian forcing uncertainty
    #assuming +/- 20% uncertainty in WMGHG forcings
    #Map the TCR cdf to the forcing scaling cdf using a beta function cdf
    beta_cdf_sf = root(lambda var: 0.05 - beta.cdf(0.5 - 1.0/3.0,var,var),x0=2.0 ).x[0]
    cdfs_gauss = 1.0 - beta.cdf(cdfs,beta_cdf_sf,beta_cdf_sf)
    sf_gauss =   (np.sum(rf_2011_mid_a[1:-2]) + np.sqrt(2.0)*erfinv(2*cdfs_gauss-1)*sig_tot)/ np.sum(rf_2011_mid_a[1:-2])
    
    return sf_gauss

# -------------------------------------------------------------------- #

# Aerosol scaing factor calculated by calculating the aerosol-only temperature response in 2017, temperature 
#     response from all forcers except aerosols in 2017, and 0.5 cdf value of warming in 2017.
#     find the required scaling factor on aerosol-only warming to make the value equal to 0.5 cdf warming in 2017 
#     minus the contribution from all other forcers (but not aerosols) in 2017.
def create_aerosol_scalingfactors(cdfs, TCR, ECS, rf_scens, sf_gauss, rf_cont, rf_nat):
    """Create the aerosol forcing scaling factors that are needed to
    ensure that warming is equal to the 0.5 cdf value in 2017 for all TCRs"""

    # Scale aerosol contribution to the correct fraction of total warming in 2017.
    # Total temperature is quoted as the values predicted by FaIR model output in 2017, for smooth join to plumes.
    # These arise because fair runs to match anthro forcing, but total temp includes natural forcing contribution.
    # Also the 4 dataset mean temperature observations we fit aren't the same as the ones used in haustein et al. 2017 nature paper. 
    temp_2017 = [0.835, 1.045, 1.260]
    # find standard deviation from temperature range
    sig_temp = (temp_2017[0] - temp_2017[1])/(erfinv((0.5-(1.0/3.0))*2 - 1)*np.sqrt(2.0))
    # Compute the TCR and ECS intervals that are consistent 
    # with likely above and likely below
    temp_2017 = sig_temp * np.sqrt(2.0) * erfinv(2*(cdfs) -1.0) + (temp_2017[1] )
    
    # prepare arrays to store output
    sf_aero = np.zeros(len(TCR),dtype=np.float)
    rf_aero_2011 = np.zeros(len(TCR),dtype=np.float)
    rf_tot_2011 = np.zeros(len(TCR),dtype=np.float)
    
    for i in range(0,len(TCR)):
        # Calculate the aerosol-only forcing response 
        aero_conc, aero_temp  = fair_scm(other_rf=rf_scens['RCP8.5']['rf_aero'],
                                         tcrecs=np.array([float(TCR[i]),float(ECS[i])]),
                                         F_2x=a*np.log(2),d=d)

        # Calculate the response to all other forcings 
        oth_forc = rf_scens['RCP8.5']['rf_co2'] + sf_gauss[i]*rf_scens['RCP8.5']['rf_gauss'] +\
                    rf_scens['RCP8.5']['rf_bc'] + rf_cont + rf_nat
        other_conc, other_temp  = fair_scm(other_rf= oth_forc,tcrecs=np.array([float(TCR[i]),float(ECS[i])]),
                                           F_2x=a*np.log(2),d=d)
        
        # rebaseline to common reference period
        aero_temp = express_as_anom(aero_temp,years) 
        other_temp = express_as_anom(other_temp,years) 
    
        # calculate aerosol scaling factor to 'fill the gap' between all other forcing and total temp. response. 
        x = (temp_2017[i] - other_temp[years==2017]) / aero_temp[years==2017]
        sf_aero[i] = np.copy(x[0])
    
        # calculate 2011 aerosol RF values and 2011 total RF values and output
        rf_aero_2011[i]= sf_aero[i] * rf_scens['RCP8.5']['rf_aero'][years==2011]
        forc = rf_scens['RCP8.5']['rf_co2'] + sf_gauss[i]*rf_scens['RCP8.5']['rf_gauss'] +\
                rf_scens['RCP8.5']['rf_bc'] + sf_aero[i]*rf_scens['RCP8.5']['rf_aero'] +\
                rf_cont
        rf_tot_2011[i] = forc[years==2011]
        
    return sf_aero, rf_aero_2011, rf_tot_2011

# -----------------------------------------------
