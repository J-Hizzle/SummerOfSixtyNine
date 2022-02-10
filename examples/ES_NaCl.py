"""
Example script to showcase the functionality of the implemented Ewald summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell
"""
# %%
import SummerOfSixtyNine.ssn as ssn
from SummerOfSixtyNine.ssn.structure import Structure
# %%
# specify input parameters
real_cut    = 1       # arbitrary for ES
outval      = 'E'       # calculate energy
technique   = 'DS'      # use direct sum method
sigma       = 1         # define Gaussian parameter
n_cut       = 1
k_cut       = 1


# specify struc parameters formula, oxidation_states, crystal_type, cell_parameters
formula             = 'NaCl'
atomic_numbers      = [11, 17]
oxidation_states    = [1, -1]
crystal_type        = 'NaCl'
cell_parameters     = [5.64e-10, 4]

# instantiate struc class
struc = Structure(formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters)
# %%
E_real = ewald_real(struc, sigma, real_cut)
print(E_real)
# %%
