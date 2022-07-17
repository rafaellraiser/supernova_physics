# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:33:07 2022

@author: Rafael de Lima Raiser
"""

#Simulate neutrino interactions with LVD and generate histograms
#Configurations
#Fiducial mass: 1 kton cintillator
#N_p=1.7e+32
#N_e=6.3e+32
#N_C=7.7e+31
#HTE (high threshold energy): 4 MeV (internal tanks)
#LTE: 0.8 MeV
#Operational time: 30 yrs

# FUNCTIONS

# Generate random energies acording to Garching flux
def random_gen_garching(par):
    import numpy as np
    import math
    F0=par[0]
    E0=par[1]
    a=par[2]
    
    f=lambda E: ( (1+a)**(1+a)*F0/(math.gamma(1+a)*E0) )* (E/E0)**a *np.exp(-(1+a)*E/E0)
    Emax=100
    xmax=E0*a/(1+a)
    M=f(xmax)
    x=Emax*np.random.rand()
    y=M*np.random.rand()
    while y>f(x):
        x=Emax*np.random.rand()
        y=M*np.random.rand()
    return x

# Model for the Flux of neutrinos - Garching model
def flux(E,par):
    import math
    import numpy as np
    F0=par[0]
    E0=par[1]
    a=par[2]
    
    f=( (1+a)**(1+a)*F0/(math.gamma(1+a)*E0) )*(E/E0)**a * np.exp(-(1+a)*E/E0)
    
    return f

# CODE

import numpy as np
import scipy.integrate as integrate
from matplotlib import pyplot as plt
import math

# Detector parameters
#N_p=1.7e+32 # number of p
#N_e=6.3e+32 # number of e-
N_C=7.7e+31 # number o 12c
T=100 *(365*24*3600 )  #operational time in s
HTE=4 # /MeV
LTE=0.8 # /MeV

#Flux parameters [total flux /cm^-2s-1, avg energy /MeV, pinching]
par_e=[10,10,2.3]
par_ebar=[10,12,2.3]
par_x=[10,15,2.3] # flux of mu+tau neutrinos (antineutrinos assumed to be the same)


# Cross sections for IBD and NC scattering on 12C
## NC: vx + C -> vx+ *C for neutrinos
from cross_sections import cs_NC_C_nu
## NC: vx + C -> vx+ *C for anti-neutrinos
from cross_sections import cs_NC_C_nubar

## IBD
#from cross_sections import cs_ibd

# Average number of interactions
## Parameters of added fluxes are estimated
import sum_garching

#flux_IBD=lambda E: flux(E,F0[1],E0[1],a[1])
#par_IBD=[F0[1],E0[1],a[1]]

flux_NCnu=lambda E: flux(E,par_e) + flux(E,par_x)
par_NCnu=sum_garching(par_e,par_x)

flux_NCnubar= lambda E: flux(E,par_ebar) + flux(E,par_x)
par_NCnubar=sum_garching(par_ebar,par_x)

## IBD
#Navg_IBD,_=integrate.quad(lambda E: N_p*T*flux_IBD(E)*cs_ibd(E),0,np.Inf)

## NC nu (nu_e,nu_x)
Navg_NC_nu,_=integrate.quad(lambda E: N_C*T* flux_NCnu(E) *cs_NC_C_nu(E),0,np.Inf)

## NC nubar (nu_ebar,nu_xbar)
Navg_NC_nubar,_=integrate.quad(lambda E: N_C*T* flux_NCnubar(E)*cs_NC_C_nubar(E),0,np.Inf)


# ITERACTION OVER RUNS
runs=100

Nev_NC_nu=np.zeros(runs)
Nev_NC_nubar=np.zeros(runs)

for run in range(runs):
    
    # Sample k events from average using poisson distribution
    
    #k_IBD=np.random.poisson(Navg_IBD)
    k_NC_nu=np.random.poisson(Navg_NC_nu)
    k_NC_nubar=np.random.poisson(Navg_NC_nubar)
        
    #Ev_IBD=np.zeros(k_IBD)
    Ev_NC_nu=np.array([])
    Ev_NC_nubar=np.array([])
    
    for i in range(k_NC_nu):
        r=random_gen_garching(par_NCnu)
        if r>15.11:
            Ev_NC_nu=np.concatenate( (Ev_NC_nu,[r]) )
    
    for i in range(k_NC_nubar):
        r=random_gen_garching(par_NCnubar)
        if r>15.11:
            Ev_NC_nubar=np.concatenate( (Ev_NC_nubar,[r]) )

    #  Number of events
    Nev_NC_nu[run]=len(Ev_NC_nu)
    Nev_NC_nubar[run]=len(Ev_NC_nubar)
    
# Average number of events and std

Nev_NC=Nev_NC_nu+Nev_NC_nubar
Nev_NC_avg=np.mean(Nev_NC)
Nev_NC_std=np.std(Nev_NC)




