# Functions describing efficiencies, thresholds and resolutions of detectors

# Super-Kamiokande

# th_*(x) returns 1 if x is inside the acceptable interval of detection
# res_*(x) returns the gaussian resolution for energy x
# cut_eff_*(x) returns the probability of a particle detected with energy
# eff_*(x): efficiency

# Super-Kamiokande

# E=np.array
def eff_diffuse_SKIII(Ee):
    
    if Ee<16:
        eff=0
    elif Ee<17:
        eff=0.55
    elif Ee<18:
        eff=0.61
    elif Ee<19:
        eff=0.72
    elif Ee<20:
        eff=0.76
    elif Ee<24:
        eff=0.82
    else:
        eff=0.92

    return eff

# E=np.array
def eff_diffuse_SKII(Ee):
    if Ee<17.5:
        eff=0
    elif Ee<18:
        eff=0.47
    elif Ee<19:
        eff=0.52
    elif Ee<20:
        eff=0.56
    elif Ee<21:
        eff=0.69
    elif Ee<26:
        eff=0.74
    else:
        eff=0.83
    return eff 

# Source: K. Abe et al arXiv:1606.07538
def res_SKIV(Ee,Ed):
    import numpy as np
    std= -0.0839+ 0.349*np.power(Ed,0.5) +0.0397 *Ed
    
    y=np.power(std*np.sqrt(2*np.pi),-1)*np.exp(-0.5*np.power(Ee-Ed,2)*np.power(sigma,-1))
    return y

# Source: K. Abe et al arXiv:1606.07538
def res_SKIII(Ee,Ed):
    import numpy as np
    std= -0.123+ 0.376*np.power(Ed,0.5) +0.0349 *Ed
    
    y=np.power(std*np.sqrt(2*np.pi),-1)*np.exp(-0.5*np.power(Ee-Ed,2)*np.power(std,-1))
    return y

# Normalized gaussian
def res_SKIII(Ee,Ed):
    import numpy as np
    std= -0.123+ 0.376*np.power(Ed,0.5) +0.0349 *Ed
    y=np.power(std*np.sqrt(2*np.pi),-1)*np.exp(-0.5*np.power(Ee-Ed,2)*np.power(std,-1))
    return y

# Ee como array

def eff_SKIII(Ee):

    import scipy.integrate as integrate
    import numpy as np
    Ee=np.linspace(0,100,1000)
    Ee_th=3.5 # threshold /MeV

    std= -0.123+ 0.376*np.power(Ee,0.5) +0.0349 *Ee
    gauss=lambda x: (std[n]*np.sqrt(2*np.pi))**(-1)*np.exp(-0.5*(Ee[n]-x)**2 / std[n]**2 )    
    error=np.zeros(len(Ee))
    eff=np.zeros(len(Ee))

    for n in range(len(Ee)):
    
        eff[n],error[n]=integrate.quad(gauss,Ee_th,np.Inf)
        
    return eff,max(error)
    