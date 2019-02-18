## Python code to import and process the different historical temperature observation datasets used in Chapter 1, SR1.5. 

# Written by Stuart Jenkins (stuart.jenkins@wadham.ox.ac.uk) (18/12/2018)

# Import statements
import numpy as np
import scipy as sp
from statsmodels.api import OLS
import statsmodels.tools.tools
from pandas import DataFrame

from Support_codes.fair_scm import *
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


# -------------------------------------------------
# Import and rebaseline the observations ready for plotting
# -------------------------------------------------
def temp_import():
    
    """
    Imports the HadCRUT4, HadCRUT4-CW, NOAA and GISTEMP datasets, re-baselines them to 1850-1900
    """
    
    # define the baseline year range, and common reference range
    base_low=1850.
    base_high=1900.
    com_ref_low=1880.
    com_ref_high=2017.
    # define variable representing the frequency of temperature observations data ('mon' = monthly)
    temp_freq='mon'

    # -------------------------------------------------
    ## Import the temperature observation datasets ##
    #Specify the GMST best-estimate temperature timeseries files to load from
    gmst_files = {'HadCRUT4':'./Data/HadCRUT.4.6.0.0.monthly_ns_avg.txt',
    'GISTEMP':'./Data/GLB.Ts+dSST.csv',
    'NOAA':'./Data/aravg.mon.land_ocean.90S.90N.v4.0.1.201803.asc',
    'Cowtan-Way':'./Data/had4_krig_v2_0_0.txt'}

    gmst_names = gmst_files.keys()
    # make a common years vector, which we can use as the years variable on all imported temperature datasets
    years_com = np.arange(1850. + 1./24,1850. + 1./24 + (2020)*1./12,1.0/12)[:-1]

    # define dictionary gmst to hold the temperature data and its averages etc.
    gmst = {}

    # Go through the datasets imported from the files referenced in 'gmst_files' above and load them
    for key in gmst_names:

        if key in ['HadCRUT4','Cowtan-Way']:
            data = np.genfromtxt(gmst_files[key])
            temps = data[:,1]
            years = years_com[:len(temps)]

        if key in ['GISTEMP']:
            f_giss = open(gmst_files[key],'r')
            temps = []
            counter = 0
            for line in f_giss:
              if counter>=2:
                  temps.extend([float(f) for f in line.split(',')[1:13] if f != '***'])
              counter = counter + 1
            temps=np.array(temps)
            years = years_com[years_com>1880.][:len(temps)]

        if key in ['NOAA']:
            data = np.genfromtxt(gmst_files[key])
            temps = data[:,2]
            years = years_com[years_com>1880.][:len(temps)]


        gmst[key] = {'Temp':temps,'Years':years}

    #Set the datasets to a common reference period        
    hc_ref = np.mean(gmst['HadCRUT4']['Temp'][np.logical_and(gmst['HadCRUT4']['Years']>=com_ref_low,
                    gmst['HadCRUT4']['Years']<(com_ref_high+1))]) - np.mean(gmst['HadCRUT4']['Temp'][np.logical_and(gmst['HadCRUT4']['Years']>=base_low,
                                                    gmst['HadCRUT4']['Years']<(base_high+1))])
    for key in gmst_names:
        gmst[key]['Temp'] = gmst[key]['Temp'][gmst[key]['Years'] < 2018.]
        gmst[key]['Years'] = gmst[key]['Years'][gmst[key]['Years'] < 2018.]
        #Express relative to a common base period
        gmst[key]['Temp'] = gmst[key]['Temp'] - np.mean(gmst[key]['Temp'][np.logical_and(gmst[key]['Years']>=com_ref_low,
                                                                  gmst[key]['Years']<(com_ref_high+1))])
        #Set NOAA and GISTEMP datasets relative to HadCRUT4 value over the base period 
        if key in ['NOAA','GISTEMP']:
            gmst[key]['Temp'] = gmst[key]['Temp'] + hc_ref
        else: 
            gmst[key]['Temp'] = gmst[key]['Temp'] - np.mean(gmst[key]['Temp'][np.logical_and(gmst[key]['Years']>=base_low,gmst[key]['Years']<(base_high+1))])

    return gmst

# -------------------------------------------------


# -----------------------------------------------
# Find the min, mean and max values from the temperautre observations
# -----------------------------------------------
def calc_mean_min_max(gmst):

    """
    Requires gmst to have dictionary strings: HadCRUT4, Cowtan-Way, GISTEMP, NOAA
    """
    
    obs_max = np.zeros_like(gmst['HadCRUT4']['Years'])
    obs_min = np.zeros_like(gmst['HadCRUT4']['Years'])
    obs_mean = np.zeros_like(gmst['HadCRUT4']['Years'])
    for y in range(0,len(gmst['HadCRUT4']['Years'])): 
        year_vals = []
        #Loop over AR5 datasets and Cowtan-Way
        for ob in ['HadCRUT4','NOAA','GISTEMP','Cowtan-Way']:
            # collect the temperature value at a given year in each dataset and store in val
            val = gmst[ob]['Temp'][gmst[ob]['Years']==gmst['HadCRUT4']['Years'][y]]
            if len(val)>0:
                year_vals.append(val)
        # find the min, mean and max values from each year
        obs_max[y] = np.max(year_vals)
        obs_min[y] = np.min(year_vals)
        obs_mean[y] = np.mean(year_vals)

    # save as entries in gmst
    gmst['Temp-max'] = obs_max
    gmst['Temp-min'] = obs_min
    gmst['Temp-mean'] = obs_mean
    
    return gmst

# -------------------------------------------------


# -----------------------------------------------
# Using OLS regression to scale anthropogenic and natural contributions to observed GMST data
# Methodology follows Haustein et al. (Scientific Reports, 2017)
# -----------------------------------------------
def calc_gwi(obs,obs_years,reg_type='mon',base_low=1850.,base_high=1900, name=''):
    
    #Express the observations relative to the base period 
    obs = obs - np.mean(obs[np.logical_and(obs_years>=base_low,obs_years<(base_high+1))])

    #Load the best estimate forcings from Piers
    forc_file = './Data/Annualforcings_Mar2014_GHGrevised.txt'
    data = np.genfromtxt(forc_file,skip_header=4)
    years = data[:,0]
    tot_forc = data[:,13]
    ant_forc = data[:,14]
    
    #Integrate anthropogenic and natural forcing with standard FAIR parameters
    C, t_nat = fair_scm(other_rf=tot_forc-ant_forc)
    C, t_anthro = fair_scm(other_rf=ant_forc)
    #Express relative to the centre of the base period
    t_nat = t_nat - np.mean(t_nat[np.logical_and(years>=base_low,years<base_high+1)])
    t_anthro = t_anthro - np.mean(t_anthro[np.logical_and(years>=base_low,years<base_high+1)])
    # -----------------------------------------------
    
    
    # Prepare the temperatures run through FaIR, so they lie on same year-grid as observations, so they can be compared
    # -----------------------------------------------
    #Interpolate the annual forced responses to the grid of the observed data
    if reg_type !='mon':
        t_nat = np.interp(obs_years+0.5, years+0.5, t_nat)
        t_anthro = np.interp(obs_years+0.5, years+0.5, t_anthro)
    else:
        t_nat = np.interp(obs_years, years+0.5, t_nat)
        t_anthro = np.interp(obs_years, years+0.5, t_anthro)

    #Linearly project the final half year
    t_anthro[obs_years>(years[-1]+0.5)] = 12*(t_anthro[obs_years<=(years[-1]+0.5)][-1] - t_anthro[obs_years<=(years[-1]+0.5)][-2]) * (obs_years[obs_years>(years[-1]+0.5)] - obs_years[obs_years<=(years[-1]+0.5)][-1]) \
    +t_anthro[obs_years<=(years[-1]+0.5)][-1]
    t_nat[obs_years>(years[-1]+0.5)] = 12*(t_nat[obs_years<=(years[-1]+0.5)][-1] - t_nat[obs_years<=(years[-1]+0.5)][-2]) * (obs_years[obs_years>(years[-1]+0.5)] - obs_years[obs_years<=(years[-1]+0.5)][-1]) \
    +t_nat[obs_years<=(years[-1]+0.5)][-1]
    # -----------------------------------------------
    
    #Use scipy defined OLS regression function to complete OLD regression of observations data on natural and anthropogenic warming with a constant
    y = np.copy(obs)
    x = DataFrame({'x1': (t_anthro), 'x2': (t_nat)})
    # add constant vector on to dataframe we will fit to temp observations
    x = statsmodels.tools.tools.add_constant(x)
    # complete OLS regression of anthropogenic and natural temperatures (found from FaIR integrated best estimate forcing) onto given observed temperature dataset.
    model = OLS(y, x)
    result = model.fit()
    # collect output scaling factors for anthro and natural temperature timeseries
    sf = result.params

    #Form scaled anthropgenic warming index
    awi = t_anthro * sf['x1']
    #Scaled natural warming index
    nwi = t_nat * sf['x2']
    #Scaled total externally forced warming index
    gwi = awi + nwi
    
    print(name, ' AWI scale factor: ', sf['x1'], '\n', name, ' NWI scale factor: ', sf['x2'])

    
    return awi, nwi

# -------------------------------------------------
