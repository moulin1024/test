import numpy as np
import matplotlib.pyplot as plt

def gaussian_wake_model(location,wt):
    # A implementation of the Gaussian wake model based on the following papers: 
    # Bastankhah, Majid, and Fernando Porté-Agel. "Experimental and theoretical study of wind turbine wakes in yawed conditions." Journal of Fluid Mechanics 806 (2016): 506-541.
    
    # a emprical fitting of wake expansion rate provided by 
    # Niayifar, Amin, and Fernando Porté-Agel. "Analytical modeling of wind farms: A new approach for power prediction." Energies 9.9 (2016): 741.
    k = (0.3837*wt.ti_inflow+0.003678)
    # Model constants provided by Bastankhah & Porté-Agel
    alpha = 2.32 
    beta = 0.154

    # Change the reference point to the given wind turbine
    x = location[0]-wt.location[0]
    y = location[1]-wt.location[1]
    z = wt.hub_height

    theta = 0.3*wt.gamma/np.cos(wt.gamma)*(1-np.sqrt(1-wt.CT*np.cos(wt.gamma)))
    # Onset of the far wake region (x0)
    x0 = wt.diameter*(np.cos( wt.gamma)*(1+np.sqrt(1- wt.CT)))/(np.sqrt(2)*(alpha*wt.ti_inflow+beta*(1-np.sqrt(1- wt.CT))))
    
    sigma_y= k*(x-x0)+wt.diameter*np.cos( wt.gamma)/np.sqrt(8)
    sigma_z= k*(x-x0)+wt.diameter/np.sqrt(8) 
    # Mask the near-wake region (x<x0) where the analytical model doesn't apply 
    if x>x0:
        # Eq. 7.4
        tmp1 = (1.6+np.sqrt( wt.CT))*(1.6*np.sqrt(8*sigma_y*sigma_z/(wt.diameter**2*np.cos( wt.gamma)))-np.sqrt( wt.CT)) 
        tmp2 = (1.6-np.sqrt( wt.CT))*(1.6*np.sqrt(8*sigma_y*sigma_z/(wt.diameter**2*np.cos( wt.gamma)))+np.sqrt( wt.CT)) 

        delta = theta*x0+wt.diameter*theta/14.7*np.sqrt(np.cos( wt.gamma)/(k**2* wt.CT))*((2.9+1.3*np.sqrt(1- wt.CT)- wt.CT)*np.log(tmp1/tmp2))
        
        exp_part = (np.exp(-0.5*((y-delta)/sigma_y)**2)*np.exp(-0.5*((z-wt.hub_height)/sigma_z)**2))
        # Eq. 7.1 
        delta_u = -wt.u_inflow*(1-np.sqrt(1-wt.CT*np.cos(wt.gamma)/(8*(sigma_y*sigma_z/wt.diameter**2))))*exp_part
    else:
        delta_u = 0

    return delta_u

def turbulent_intensity_model(location,wt):
    # An implementation of the bi-gaussian emprical model by 
    # Qian, G. W., & Ishihara, T. (2018). A New Analytical Wake Model for Yawed Wind Turbines. Energies, 11(3), 665.
    gamma = wt.gamma
    CT_BP = wt.CT*np.cos(gamma)**2
    wt_location = wt.location
    zh = wt.hub_height
    d = wt.diameter
    I_inf = wt.ti_inflow
    # Qian and Ishihara adopted a different definition for CT
    CT = CT_BP/np.cos(gamma)**2 
    CTp = CT*np.cos(gamma)**3 
    theta0 = 0.3*gamma/np.cos(gamma)*(1-np.sqrt(1-CTp)) 
    
    x = location[0]-wt_location[0]
    y = location[1]-wt_location[1]
    z = zh
    
    kstar = 0.11*CTp**1.07*I_inf**0.2 
    epsilonstar = 0.23*CTp**(-0.25)*I_inf**0.17 

    sigma = kstar*x+epsilonstar*d 
    
    if x>=0:
        # Exception for zero yaw angle in which the original model will divide zero 
        if  gamma-0>1e-6:
            sigma0 = np.sqrt((CT*np.cos(gamma)**2*(np.sin(gamma)+1.88*np.cos(gamma)*theta0))/(44.4*theta0))*d 
            x0 = (sigma0 - epsilonstar*d)/kstar 

            if x<=x0:
                delta = theta0*x 
            else:
                temp1 = (sigma0/d+0.2*np.sqrt(CTp))*(sigma/d-0.2*np.sqrt(CTp)) 
                temp2 = (sigma0/d-0.2*np.sqrt(CTp))*(sigma/d+0.2*np.sqrt(CTp)) 
                delta = theta0*x0+d*(np.sqrt(CT*np.cos(gamma))*np.sin(gamma))/(18.24*kstar)*np.log(abs((temp1)/(temp2)))
        else:
            delta = 0

        r = np.sqrt((y-delta)**2+(z-zh)**2) 


        dd = 2.3*CTp**(-1.2) 
        e = 1.0*I_inf**(0.1) 
        q = 0.7*CTp**(-3.2)*I_inf**(-0.45)/(1+x/d)**2 

        if r/d<=0.5:
            k1 = np.cos(np.pi/2*(r/d-0.5))**2 
            k2 = np.cos(np.pi/2*(r/d+0.5))**2 
        else:
            k1 = 1 
            k2 = 0

        phi2 = k1*np.exp(-0.5*((r-d/2)/sigma)**2) + k2*np.exp(-0.5*((r+d/2)/sigma)**2) 

        G = 1/(dd+e*x/d+q) 
        delta_I = G*phi2 

    if x<0:
        delta_I = 0
    I_new = delta_I 
    return I_new