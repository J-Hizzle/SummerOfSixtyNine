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
import csv
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
# set file to write generated data to
csvfile = '../data/ewald_sum_results.csv'
# %%
# open the file in the write mode
f = open(csvfile, 'w')

# create the csv writer
writer = csv.writer(f)
writer.writerow(['n_cut', 'k_cut', 'sigma', 'energy', 'madelung', 'time', 'memory'])
f.close()# make array of energies depending on cutoff values
# %%
n_cut_vals = np.arange(5, 10)
k_cut_vals = np.arange(5, 25)
eta = 8 * 1/np.sqrt(2) * 1/(structure.volume**2)**(1/3) * np.pi
#sig = 1/eta
#sig = 1e-10
#print(sig)
sigma_vals = np.linspace(1e-11, 1e-9)
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
            
            # open csv file in append mode
            f = open(csvfile, 'a')
            # create the csv writer
            writer = csv.writer(f)

            # append values to csv file
            writer.writerow([n_cut, k_cut, sigma, E, M, time_total, mem_use])

            # close the file
            f.close()

        
        timings_n.append(timings_k)
        mem_uses_n.append(mem_uses_k)
        energies_n.append(energies_k)
        madelungs_n.append(madelungs_k)
    
    timings.append(timings_n)
    mem_uses.append(mem_uses_n)
    energies.append(energies_n)
    madelungs.append(madelungs_n)