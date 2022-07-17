# MSW mixing of flavors in supernova
# Adiabatic transitions at low density and high density layers
# f are functions of energy
def mix_ad(f_e,f_ebar,f_x,component='nu_ebar',hierarchy='nh'):
    import numpy as np
    p_nh=0.0212
    pbar_nh=0.6783
    p_ih=0.3005
    pbar_ih=0.0212
    
    if hierarchy=='nh':
    
        if component=='nu_ebar':
            f=lambda E: pbar_nh*f_ebar(E)+(1-pbar_nh)*f_x(E)
        elif component=='nu_e':
            f=lambda E: p_nh*f_e(E)+(1-p_nh)*f_x(E)
        elif component=='nu_x':
            f=lambda E: 0.25*( (1-p_nh)*f_e(E)+ (1-pbar_nh)*f_ebar(E) + (2+p_nh+pbar_nh)*f_x(E)  )
    
    elif hierarchy=='ih':
        if component=='nu_ebar':
            f=lambda E: pbar_ih*f_ebar(E)+(1-pbar_ih)*f_x(E)
        elif component=='nu_e':
            f=lambda E: p_ih*f_e(E)+(1-p_ih)*f_x(E)
        elif component=='nu_x':
            f=lambda E: 0.25*( (1-p_ih)*f_e(E)+ (1-pbar_ih)*f_ebar(E) + (2+p_ih+pbar_ih)*f_x(E)  )
   
    return f
    
        