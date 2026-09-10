import numpy as np # numpy provides powerful multi-dimensional arrays to hold and manipulate data
import matplotlib.pyplot as plt # matplotlib provides powerful functions for plotting figures
import matplotlib 
import astropy.units as u # astropy provides unit system and constants for astronomical calculations
import astropy.constants as const
from IPython.display import Latex # import Latex module so we can display the results with symbols
 

def Read(filename):
    """ This function will be used to read the files infomation line by line so its easier to use in future projects
        Inputs: filename - the name of the file so it can be open and read for the Read function to be used
        Outputs: time, number ,data 
    """
    file = open(filename,'r') # this will open and read the file 
    
    line1 = file.readline() # read the first line 
    label,value = line1.split() # will split each label and value
    time = float(value)*u.Myr # converting units for time
    line2 = file.readline() # read the second line 
    label,value = line2.split() # again just splitting the labels and values
    number = float(value) # the number will be a float from what we calculated from the value function
# close file 
    file.close()
#separates the data types and skips the first 3 lines
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    return time, number, data # this will just return the time number and data when it the reader asks to use this function in projects 