import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants
from ReadFile import Read
 # this is so we can read the file and take the data for our calculations

def ComponentMass(filename, particle_type):
    """ This function will return the total mass of any galaxy component the reader decides. We will be using the read file from the last assignment to help us use the read fucnction.
        The mass will be returned in units of 10*12 solar mass while also rounding up to 3 decimal places.
    Inputs: filename - the name of the file so it can be open and read for the Read function to be used
        particle type is the type of particle (Type 1 = Dark Matter, Type 2 = Disk Stars, Type 3 = Bulge Stars)
    Outputs: Totmass (solar mass)- total mass rounded 3 places and in units of 10^12 solar mass
    """
    #read the data from the imported file
    time, number,data = Read(filename)
    # this is used to go where the particle type is located in the data easier (1-3)
    types = data['type'] == particle_type
    i = np.where(types)[0] # making sure it goes to the location for the info
    # this is used to calculate the total mass of the galaxy component
    mass = np.sum(data['m'][i])  
    # we are going to round the mass before it finally return the total mass
    Totmass = np.round(mass* 1e-2,3)  
    return Totmass

table = [] # this is the start of the table so it can be easily calculated
print("Galaxy Name | Halo Mass (Mo) | Disk Mass (Mo) | Bulge Mass (M0) | Total (Mo) | fbar|") # this will print the table so easier to locate each item
filename = "M33_000.txt" # the selected file for this table
halo = ComponentMass(filename,1)  # halo stars calculated
table.append([halo]) # appends the halo calculations
disk = ComponentMass(filename,2) # the disk stars to be calculated
table.append([disk])  # append the disk calcualtions
bulge = ComponentMass(filename,3) # the bulge calculations
table.append([bulge]) # the bulge calculations being appended to the table
LocalGMass = halo + disk + bulge # this is to sum up the local group mass
table.append([LocalGMass]) # append the mass to the table
fbar = np.round(disk / LocalGMass, 3) # we need calculate the baryon fraction is the total stellar mass over total mass
table.append([fbar])  # now append this to the table
print(table) # print the table 