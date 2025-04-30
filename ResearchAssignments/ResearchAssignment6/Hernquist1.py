import numpy as np # numpy provides powerful multi-dimensional arrays to hold and manipulate data
import matplotlib.pyplot as plt # matplotlib provides powerful functions for plotting figures
import matplotlib 
import astropy.units as u # astropy provides unit system and constants for astronomical calculations
import astropy.constants as const

def Hernquist(r, Mhalo, a):
    """ This function will compute the density for radii within, using the theoretical profile.  
        Inputs: r: `(kpc) this is the distance from the center  
        a: (kpc) the radius length  of the mass density's movement
        Mhalo: (solar mass) the total mass of the halo 
        Returns: HernDens:  the density distribution 
    """ 
    x = ((Mhalo*u.Msun) /(2*np.pi)) # first part of the equation and make conversions
    y = ((a*u.kpc) / ((r*u.kpc)*(r*u.kpc + a*u.kpc )**3)) # the second part of the equation
    HernDens = x*y # we now multiple both equations 
    return HernDens  # this is the hernquist density profile in solar mass units