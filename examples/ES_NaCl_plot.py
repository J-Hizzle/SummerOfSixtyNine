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
from scipy.constants import N_A


import numpy as np
import SummerOfSixtyNine.ssn as ssn
from SummerOfSixtyNine.ssn.structure import Structure
# %%
# specify input parameters
outval      = ['E', 'M']       # calculate E and M
technique   = 'ES'      # use direct sum method

# specify structure parameters formula, oxidation_states, crystal_type, cell_parameters
formula             = 'NaCl'
atomic_numbers      = [11, 17]
oxidation_states    = [1, -1]
crystal_type        = 'NaCl'
cell_parameters     = [5.64e-10, 4]

# instantiate structure class
structure = Structure(formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters)
# %%
# make array of energies depending on cutoff values
n_cut_vals = np.arange(20, 21)
k_cut_vals = np.arange(30, 31)
eta = 8 * 1/np.sqrt(2) * 1/(structure.volume**2)**(1/3) * np.pi
#sig = 1/eta
#sig = 1e-10
#print(sig)
sigma_vals = np.linspace(5e-8, 1e-6)
# %%
energies = []
timings = []
mem_uses = []
madelungs = []
for sigma in sigma_vals:
    timings_n = []
    mem_uses_n = []
    energies_n = []
    madelungs_n = []

    for n_cut in n_cut_vals:
        timings_k = []
        mem_uses_k = []
        energies_k = []
        madelungs_k = []

        for k_cut in k_cut_vals:
            # measure time before running ssn
            start_time = time.time()
            
            # run ssn while measuring memory usage
            mem_return = mp.memory_usage((ssn.run_ssn, [structure, n_cut, technique, outval, k_cut, sigma]), interval=.001, retval=True)
            
            # measure time after running ssn
            time_total = (time.time() - start_time)

            # get memory usage and E from memory_usage return values
            mem_use = np.max(mem_return[0]) - np.min(mem_return[0])
            E = mem_return[1][0]
            M = mem_return[1][1]

            print('M({0},{1},{2}) ='.format(sigma, n_cut, k_cut), M)

            # append values to lists
            timings_k.append(time_total)
            mem_uses_k.append(mem_use)
            energies_k.append(E)
            madelungs_k.append(M)
        
        timings_n.append(timings_k)
        mem_uses_n.append(mem_uses_k)
        energies_n.append(energies_k)
        madelungs_n.append(madelungs_k)
    
    timings.append(timings_n)
    mem_uses.append(mem_uses_n)
    energies.append(energies_n)
    madelungs.append(madelungs_n)
 # %%
timings = np.array(timings)
mem_uses = np.array(mem_uses)
energies = np.array(energies)
madelungs = np.array(madelungs)
# %%
print(np.shape(energies))
# %%
plt.plot(k_cut_vals[1:], madelungs[0][0][1:], 'o')
plt.xlabel(r'real space cutoff shell')
plt.ylabel(r'total electrostatic E in $10^{-6}$ kJ/mol')
plt.title(r'Divergence of the direct sum method')
plt.savefig(fname='../data/direct_sum_divergence_1', dpi=900)
# %%