"""
Example script to showcase the functionality of the implemented direct summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell
"""
# %%
import SummerOfSixtyNine.ssn as ssn
from SummerOfSixtyNine.ssn.structure import Structure
# %%
# specify input parameters
n_cut       = 10         # define real space cutoff shell
outval      = 'E'       # calculate energy
technique   = 'DS'      # use direct sum

# specify struc parameters formula, oxidation_states, crystal_type, cell_parameters
formula             = 'NaCl'
atomic_numbers      = [11, 17]
oxidation_states    = [1, -1]
crystal_type        = 'NaCl'
cell_parameters     = [5.64e-10, 4]

# instantiate struc class
structure = Structure(formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters)
# %%
E = ssn.run_ssn(structure, n_cut, technique, outval)[0]
# %%

