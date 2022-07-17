# Models describing v emission in supernovae
## Given parameters, returns the function describing spectrum

#Fermi-Dirac spectrum: F(E)=(120/7pi^4)(Etot/T^4)*E^2/(1+exp(E/T))
# free parameters: Etot: total energy, T: temperature

#Garching spectrum: F(E)= (1+a)^(1+a)Etot/(E0^2*gamma(1+a) )(E/E0)^a exp(-(1+a)E/E0)
# free parameters: E0: mean energy, Etot: total energy, a: pinching

# A FD spectrum with temperature T is ~equivalent to Garching with E0=3.15T, a=2.3

def spectrum_fd(E,Etot,T):
    import numpy as np
    f=(120/(7*np.pi**4))*(Etot/(T**4))*np.power(E,2)*np.power((1+np.exp(E/T)),-1)
    return f

def spectrum_garching(E,Etot,E0,a):
    import numpy as np
    import math
    f=(1+a)**(1+a)*Etot/(E0**2*math.gamma(1+a) )*np.power(E/E0,a) * np.exp(-(1+a)*E/E0)
    return f