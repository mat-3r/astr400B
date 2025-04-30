import numpy as np # numpy provides powerful multi-dimensional arrays to hold and manipulate data
import matplotlib.pyplot as plt # matplotlib provides powerful functions for plotting figures
import matplotlib 
import astropy.units as u # astropy provides unit system and constants for astronomical calculations
import astropy.constants as const

def NFW(r, rho, a2):
    """ This function will compute the spherically averaged density of Dark Matter halos 
        Inputs: r: `(kpc) this is the distance from the center of the Dark Matter halo to the measuring density  
        a: (kpc) the scale radius 
        rho: (solar mass) inital density 
        Returns: NFWDens:  the density distribution 
    """ 
    x1 = (rho*u.Msun/u.kpc**3) # first part of the equation and make conversions
    y2 = (r*u.kpc)/(a2*u.kpc) # the second part of the equation
    z = ((1 + (r*u.kpc)/(a2*u.kpc))**2) # the addition to the second part
    x3 = (y2) * z # take the two bottom parts together
    NFWDens = x1/x3 # we now divide both equations 
    return NFWDens  # this is the Navarro-Frank density profile in solar mass units