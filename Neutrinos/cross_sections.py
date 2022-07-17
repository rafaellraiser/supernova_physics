# Cross section for common neutrino interactions


#  Inverse beta-decay
## Eth=1.8 MeV
def cs_ibd(Ev):
    import numpy as np
    Mnp=1.293 # neutron-proton mass difference/MeV
    me=0.511  # electron mass/MeV
    Eth=Mnp+me
    d_wm= -0.00325 * (Ev-Mnp/2) #magneto-weak correction
    if Ev>Eth:
        y= 9.4e-44 * (1+d_wm)*(Ev-Mnp)*np.sqrt(np.power(Ev-Mnp,2) -me**2) # cross-sec in cm**-2
    else:
        y=0
    return y
    

# v-e scattering
## flavor={nu_e,nu_e_bar,nu_x,nu_x_bar}
def cs_ve(Ev,flavor='nu_e_bar'):
    import numpy as np
    
    G=1.166378e-11 # Fermi constant /MeV**-2
    me=0.511 # electron mass /MeV
    f=3.89e-22 # MeV^-2 to cm^2
    
    CA=-1/2
    CV=-1/2+2*0.231
    
    if flavor=='nu_e':
        A=(CV+CA+2)**2
        B=(CV-CA)**2
        C=(CV+1)**2-(CA+1)**2
    elif flavor=='nu_e_bar':
        A=(CV-CA)**2
        B=(CV+CA+2)**2
        C=(CV+1)**2 - (CA+1)**2
    elif flavor=='nu_x':
        A=(CV+CA)**2
        B=(CV-CA)**2
        C=CV**2-CA**2
    elif flavor=='nu_x_bar':
        A=(CV-CA)**2
        B=(CV+CA)**2
        C=CV**2-CA**2
 
    if Ev==0:
        cs=0
    else:
        cs=f*(G**2*me*Ev/(2*np.pi))*(A + B/3 -0.5 * C * me * np.power(Ev,-1) )
    
    return cs

# ve-16O scattering
## Eth=15 MeV
def cs_vO(Ev):
    import numpy as np
    Eth=15
    if Ev>Eth:
        z=4.7e-40 *np.power( np.power(Ev,1/4)-15**(1/4) ,6)
    else:
        z=0
    return z

# vx + 12C -> vx + 12C* (NC)
def cs_NC_C_nu(E):
    import numpy as np
    f1=[1.32236929e-50, -3.23062137e-48,  1.64010762e-46,  9.63776958e-45,-3.67035000e-43,  2.94763143e-42]
    cs=f1[0]*E**5+f1[1]*E**4+f1[2]*E**3+f1[3]*E**2+f1[4]*E+f1[5] if E>15.11 else 0
    return cs

def cs_NC_C_nubar(E):
    import numpy as np
    f1=[ 3.53573850e-51, -1.81608474e-49, -1.54443089e-46,  2.11687796e-44,-5.48641426e-43,  4.02691204e-42]
    cs=f1[0]*E**5+f1[1]*E**4+f1[2]*E**3+f1[3]*E**2+f1[4]*E+f1[5] if E>15.11 else 0
    return cs