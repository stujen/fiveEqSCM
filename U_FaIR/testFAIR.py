# ------------------------------------------------------------------
# oxfair code to test against

import numpy as np
import pandas as pd
import scipy as sp
import csv
from scipy.interpolate import interp1d
import pickle
from scipy.ndimage.filters import gaussian_filter1d as smooth
import string
from scipy.stats import linregress
import math

# define functions for within fair model

def step_conc_test(R,alpha,E,a,tau,pre_ind_C,emis2conc):
    
    E = E[:,np.newaxis]
    alpha = alpha[:,np.newaxis]
    
    R = E * emis2conc[:,np.newaxis] * a * alpha * tau * ( 1. - np.exp( -1./(alpha*tau) ) ) + R * np.exp( -1./(alpha * tau) )

    C = pre_ind_C + np.sum(R,axis=1)
    
    G_A = (C - pre_ind_C) / emis2conc
    
    return C,R,G_A

def step_forc_test(C,pre_ind_C,F_ext,f):
    
    F = np.sum(f[...,0]*np.log(C/pre_ind_C) + f[...,1]*(C - pre_ind_C) + f[...,2] * (np.sqrt(C) - np.sqrt(pre_ind_C))) + F_ext
    
    return F

def step_temp_test(S,F,q,d):
    
    S = q*F*(1-np.exp(-1/d)) + S*np.exp(-1/d)
    
    T = np.sum(S)
    
    return S,T

def g_1_test(a,tau,h):
    
    g1 = np.sum( a*tau*( 1. - (1.+h/tau)*np.exp(-100./tau) ), axis=1 )
    
    return g1

def g_0_test(a,tau,h):
    
    g0 = ( np.sinh( np.sum( a * tau * (1. - np.exp(-h/tau)) , axis=1) / g_1(a,tau,h) ) )**(-1.)
    
    return g0

def alpha_val_test(G,G_A,T,tau,a,r,h,pre_ind_C,iirf100_max = 97.0):
    
    iirf100_val = r[...,0] + r[...,1]*(G-G_A) + r[...,2]*T + r[...,3]*G_A
    
    iirf100_val = (iirf100_val>97)*97+iirf100_val*(iirf100_val<97)
        
    alpha_val = g_0_test(a,tau,h) * np.sinh(iirf100_val / g_1_test(a,tau,h))
    
    return alpha_val

def k_q_test(d,q,tcr,ecs,F_2x):
    
    k = 1.0 - (d/70.0)*(1.0 - np.exp(-70.0/d))
    
    if (tcr and ecs):
        
        q =  (1.0 / F_2x) * (1.0/(k[0]-k[1])) * np.array([tcr-k[1]*ecs,k[0]*tcr-ecs])
        
    return k, q

# define main model code

def oxfair_test(emissions,
    emis2conc = 1/(5.148*10**18 / 1e18 * np.array([12.,16.,28.]) / 28.97),
    a = np.array([[0.2173,0.2240,0.2824,0.2763],[1,0,0.,0.],[1,0,0.,0.]]),
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
    F_2x = 3.74):

    #k , q = k_q_test(d=d,q=q,tcr=1.6,ecs=2.75,F_2x=F_2x)

    G = np.cumsum(emissions,axis=1)
    C = np.zeros(emissions.shape)
    RF = np.zeros(emissions[0].shape)
    T = np.zeros(emissions[0].shape)
    alpha = np.zeros(emissions.shape)

    alpha[...,0] = alpha_val_test(G=0,G_A=0,T=0,tau=tau,a=a,r=r,h=100.,pre_ind_C=PI_C,iirf100_max = 97.0)
    C[...,0],R,G_A = step_conc_test(R = np.zeros((3,4)),alpha=alpha[...,0],E=emissions[...,0],a=a,tau=tau,pre_ind_C=PI_C,emis2conc=emis2conc)
    RF[0] = step_forc_test(C=C[...,0],pre_ind_C=PI_C,F_ext=0.,f=f)
    S,T[0] = step_temp_test(S=np.zeros(2),F=RF[0],q=q,d=d)

    for t in np.arange(1,emissions[0].size):

        alpha[...,t] = alpha_val_test(G=G[...,t-1],G_A=G_A,T=T[t-1],tau=tau,a=a,r=r,h=100.,pre_ind_C=PI_C,iirf100_max = 97.0)
        C[...,t],R,G_A = step_conc_test(R = R,alpha=alpha[...,t],E=emissions[...,t],a=a,tau=tau,pre_ind_C=PI_C,emis2conc=emis2conc)
        RF[t] = step_forc_test(C=C[...,t],pre_ind_C=PI_C,F_ext=0.,f=f)
        S,T[t] = step_temp_test(S=S,F=RF[t],q=q,d=d)
        
    return C,RF,T