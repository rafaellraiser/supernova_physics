# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 20:25:26 2022

@author: Rafael Raiser
"""

# Returns the parameters of the sum of two garching spectra
## phi: total flux
## E0: avg energy
## pinching parameter

def sum_garching(par1,par2):
    
    phi1=par1[0]
    E01=par1[1]
    a1=par1[2]
    phi2=par2[0]
    E02=par2[1]
    a2=par2[2]
    
    phi=phi1+phi2
    E0=(E01*phi1+E02*phi2)/phi
    
    aux= ( ((2+a1)/(1+a1))*E01**2*phi1 + ((2+a2)/(1+a2))*E02**2*phi2 ) / (phi*E0**2)
    
    a=(2-aux)/(aux-1)
    
    return [phi,E0,a]
