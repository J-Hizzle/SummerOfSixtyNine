"""
Example script to measure the divergence of the direct summation technique with increasing
real space cutoff values and the time and memory usage of the calculation. The results are
plotted. 
"""

# %%
from turtle import color
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import time, csv
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
csvfile = '../data/direct_sum_results.csv'
# open the file in the write mode
f = open(csvfile, 'w')

# create the csv writer
writer = csv.writer(f)
writer.writerow(['n_cut', 'energy', 'madelung', 'time', 'memory'])
f.close()# make array of energies depending on cutoff values
# %%
n_cut_vals = np.arange(0, 25)
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

     # open csv file in append mode
    f = open(csvfile, 'a')
    # create the csv writer
    writer = csv.writer(f)

    # append values to csv file
    writer.writerow([n_cut, E, M, time_total, mem_use])

    # close the file
    f.close()