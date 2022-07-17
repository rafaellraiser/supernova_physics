# Defines gravitational collapse rate function
# g(z)=c/H0 * Rsn(z) /sqrt(W_M(1+z)^3+W_l)

#Hasan Yuksel et al. “Revealing the High-Redshift Star Formation Rate with Gamma-Ray Bursts”. The Astrophysical Journal 683.1 (jul. de 2008), pp. L5–L8. ISSN: 1538-4357. DOI: 10.1086/591449. URL: http://dx.doi.org/10.1086/591449.

# inp='integrals' returns J1,J2,J3, the integration factors to correct the parameters of the diffuse flux
# inp='rate' returns the snr calculated over the array x
# inp='integrand' returns the integrand containing all terms related to the integration in redshift calculated over the array x
# rate='avg','min','max' chooses the avg, min or max value for the local snr
def snr_yuksel(x,inp,rate='avg'):
    import numpy as np
    from scipy import integrate
    # Set of cosmological parameters
    c=299792458e-3 #speed of light /km/s       
    H0=69.8 #Hubble's constant /km/s/Mpc       
    c_H0=c/H0 # Hubble radius /Mpc
    W_M=0.308 #Matter density PDG 2018
    W_L=0.692 #Dark energy density PDG 2018
    f_flux=((3.086e+24)**2 *365*24*3600)**(-1) # f_flux= Mpc^-2 yr^-1 to cm^-2 s^-1
    r0=(1.25)*1e-4 #taxa local /Mpc^-3 yr^-1
    r0min=(1.25-0.5)*1e-4
    r0max=(1.25+0.5)*1e-4
    p1=3.4
    p2=-0.3
    p3=-3.5
    k=-10
    # Chooses the avg,min or max value of the local SNR
    if rate=='avg':
        R0=r0
    elif rate=='min':
        R0=r0min
    elif rate=='max':
        R0=r0max
    
    if (inp=='integrals') and (x==0):
        
        def fn(z,n):
            import numpy as np
            y=f_flux*c_H0* R0*((1+z)**(p1*k)+((1+z)/5000)**(p2*k)+((1+z)/9)**(p3*k))**(1/k) /((1+z)**n * np.sqrt(W_M*(1+z)**3+W_L) )
            return y
        
        J1,err1=integrate.quad(lambda z: fn(z,1),0,np.Inf)
        J2,err2=integrate.quad(lambda z: fn(z,2),0,np.Inf)
        J3,err3=integrate.quad(lambda z: fn(z,3),0,np.Inf)
        
        return J1,J2,J3
    
    elif (inp=='integrand'):
        g=f_flux*c_H0* R0*np.power(np.power(1+x,p1*k)+np.power((1+x)/5000,p2*k)+np.power((1+x)/9,p3*k),(1/k) ) / np.sqrt(W_M*np.power((1+x),3)+W_L) 
        return g 
    elif (inp=='rate'):
        R=R0*np.power(np.power(1+x,p1*k)+np.power((1+x)/5000,p2*k)+np.power((1+x)/9,p3*k),(1/k) )
        return R
        
        
        