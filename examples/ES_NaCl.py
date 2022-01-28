"""
Example script to showcase the functionality of the implemented Ewald summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell
"""
# %%
from os import system
import SummerOfSixtyNine.ssn as ssn
# %%
# specify input parameters
outval      = 'E'       # calculate energy
technique   = 'ES'      # use Ewald summation
crystal     = 'NaCl'    # for sodium chloride unit cell




ssn.main()  
# %%
