"""
Example script to measure the divergence of the direct summation technique with increasing
real space cutoff values and the time and memory usage of the calculation. The results are
plotted. 
"""

# %%
from turtle import color
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
technique   = 'DS'      # use direct sum method

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
n_cut_vals = np.arange(1, 25)
# %%
energies = []
timings = []
mem_uses = []
madelungs = []

for n_cut in n_cut_vals:
    # measure time before running ssn
    start_time = time.time()
    
    # run ssn while measuring memory usage
    mem_return = mp.memory_usage((ssn.run_ssn, [structure, n_cut, technique, outval]), interval=.001, retval=True)
    
    # measure time after running ssn
    time_total = (time.time() - start_time)

    # get memory usage and E from memory_usage return values
    mem_use = np.max(mem_return[0]) - np.min(mem_return[0])
    E = mem_return[1][0]
    M = mem_return[1][1]

    print('M({0}) ='.format(np.where(n_cut_vals == n_cut)[0][0]), M)

    # append values to lists
    timings.append(time_total)
    mem_uses.append(mem_use)
    energies.append(E)
    madelungs.append(M)

# %%
timings = np.array(timings)
mem_uses = np.array(mem_uses)
energies = np.array(energies)
madelungs = np.array(madelungs)
# %%
energies_kJmol = energies * N_A * 1e-3
print(energies_kJmol)
# %%
fig, ax = plt.subplots(4, 1, figsize=(5, 20))

plt.sca(ax[0])
plt.plot(n_cut_vals[1:], energies_kJmol[1:], 'bo')
plt.ylabel(r'$E_{direct}$ in $\mathrm{kJ/mol}$')
plt.yticks((np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5)), np.around(np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5), 2))
plt.title(r'Direct sum electrostatic energies')

plt.sca(ax[1])
plt.plot(n_cut_vals[1:], madelungs[1:] , 'ro')
plt.yticks((np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5)), np.around(np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5), 5))
plt.ylabel(r'$M_{direct}$')
plt.title(r'Madelung constants')

plt.sca(ax[2])
plt.plot(n_cut_vals, timings, 'go')
plt.ylabel(r'$t_{run}$ in $\mathrm{s}$')
plt.title(r'Total runtimes')

plt.sca(ax[3])
plt.plot(n_cut_vals, mem_uses, 'o', color='purple')
plt.ylabel(r'$m$ in $\mathrm{MiB}$')
plt.title(r'Memory usages')
plt.xlabel(r'$n_{cut}$')

fig.savefig(fname='../data/Direct sum results', dpi=900)
# %%
