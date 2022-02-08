"""
Example script to showcase the functionality of the implemented direct summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell over a range of real space cutoff values with 
the results plotted in a r_cut vs. E_tot graph
"""

# %%
from matplotlib import markers
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import numpy as np
import SummerOfSixtyNine.ssn as ssn
from SummerOfSixtyNine.ssn.structure import Structure
# %%
# specify input parameters
real_cut    = 2.81e-10         # define real space cutoff value
outval      = 'E'       # calculate energy
technique   = 'DS'      # use direct sum method

# specify struc parameters formula, oxidation_states, crystal_type, cell_parameters
formula             = 'NaCl'
atomic_numbers      = [11, 17]
oxidation_states    = [1, -1]
crystal_type        = 'NaCl'
cell_parameters     = [5.64e-10, 4]

# instantiate struc class
struc = Structure(formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters)
# %%
# make array of energies depending on cutoff values
r_cut_vals = np.linspace(2.81e-10, 1e-7, 100)
# %%
energies = []
for r_cut in r_cut_vals:
    energy = ssn.run_ssn(struc, r_cut, technique, outval)
    print('E({0}) ='.format(np.where(r_cut_vals == r_cut)[0][0]), energy)
    energies.append(energy)
# %%
energies = np.array(energies)
# %%
print(energies)
# %%
plt.plot(r_cut_vals[0:26] * 1e10, energies[0:26] * 1e-24, 'o')
plt.xlabel(r'real space cutoff in $\mathrm{\AA}$')
plt.ylabel(r'total electrostatic energy in $10^{24} \mathrm{J}$')
plt.title(r'Divergence of the direct sum method')
plt.savefig(fname='../data/direct_sum_divergence', dpi=900)
# %%