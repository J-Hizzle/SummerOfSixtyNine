"""
Example script to measure the divergence of the direct summation technique with increasing
real space cutoff values and the time and memory usage of the calculation. The results are
plotted. 
"""

# %%
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import time
import memory_profiler as mp


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
r_cut_vals = np.linspace(2.81e-10, 1e-8, 10)
# %%
energies = []
timings = []
mem_uses = []
for r_cut in r_cut_vals:
    # measure time before running ssn
    start_time = time.time()
    
    # run ssn while measuring memory usage
    mem_return = mp.memory_usage((ssn.run_ssn, [struc, r_cut, technique, outval]), interval=.001, retval=True)
    
    # measure time after running ssn
    time_total = (time.time() - start_time)

    # get memory usage and energy from memory_usage return values
    mem_use = np.max(mem_return[0]) - np.min(mem_return[0])
    energy = mem_return[1]

    print('E({0}) ='.format(np.where(r_cut_vals == r_cut)[0][0]), energy)

    # append values to lists
    timings.append(time_total)
    mem_uses.append(mem_use)
    energies.append(energy)

# %%
timings = np.array(timings)
mem_uses = np.array(mem_uses)
energies = np.array(energies)
# %%
print(energies)
# %%
plt.plot(r_cut_vals * 1e10, energies * 1e15, 'o')
plt.xlabel(r'real space cutoff in $\mathrm{\AA}$')
plt.ylabel(r'total electrostatic energy in $10^{-15} \mathrm{J}$')
plt.title(r'Divergence of the direct sum method')
plt.savefig(fname='../data/direct_sum_divergence_1', dpi=900)
# %%