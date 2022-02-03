"""
Example script to showcase the functionality of the implemented direct summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell
"""

# %%
import SummerOfSixtyNine.ssn as ssn
from SummerOfSixtyNine.ssn.structure import Structure
# %%
# specify input parameters
cutoff      = 0         # only consider supercell
outval      = 'E'       # calculate energy
technique   = 'DS'      # use Ewald summation

# specify struc parameters formula, oxidation_states, crystal_type, cell_parameters
formula             = 'NaCl'
atomic_numbers      = [11, 17]
oxidation_states    = [1, -1]
crystal_type        = 'NaCl'
cell_parameters     = [5.64e-10, 4]

# instantiate struc class
struc = Structure(formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters)
# %%
energy = ssn.run_ssn(struc, cutoff, technique, outval)
# %%
