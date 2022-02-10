"""
Example script to showcase the functionality of the implemented Ewald summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell
"""
# %%
import SummerOfSixtyNine.ssn as ssn
from SummerOfSixtyNine.ssn.structure import Structure
# %%
# specify input parameters
outval      = 'E'       # calculate energy
technique   = 'ES'      # ewald summation method
sigma       = 1         # Gaussian parameter
n_cut       = 10         # Cutoff shell in real space
k_cut       = 10         # Cutoff shell in reciprocal space


# specify struc parameters formula, oxidation_states, crystal_type, cell_parameters
formula             = 'NaCl'
atomic_numbers      = [11, 17]
oxidation_states    = [1, -1]
crystal_type        = 'NaCl'
cell_parameters     = [5.64e-10, 4]

# instantiate struc class
structure = Structure(formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters)
# %%
E = ssn.run_ssn(structure, n_cut, k_cut, sigma, technique, outval)
# %%
