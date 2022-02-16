# %%
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import time
import csv
import pandas as pd
import numpy as np
from scipy.constants import N_A
# %%
csvfile = '../data/ewald_sum_results.csv'

frame = pd.read_csv(csvfile, delimiter=',')

array = frame.to_numpy()



# get indices for respective sigma and n_cut values
n_cut = 9

sigma_vals = array[:, 2][(array[:, 0] == n_cut) & (array[:, 1] == 5)]

plt.ylabel(r'$M_{ewald}$')
plt.title(r'Ewald sum convergence for different $\sigma$ ($n_{cut} = 9$))')
plt.xlabel(r'$k_{cut}$')
plt.xticks(np.arange(5, 25, 5))

for sigma in sigma_vals[1:]:
    indices = (array[:, 0] == n_cut) & (array[:, 2] == sigma)

    #n_cut_vals = array[:, 0]
    k_cut_vals = array[:, 1][indices]
    #energies = array[:, 3][indices]
    madelungs = array[:, 4][indices]
    #timings = array[:, 5][indices]
    #mem_uses = array[:, 6][indices]

    #n_cut_vals = np.array(n_cut_vals)
    k_cut_vals = np.array(k_cut_vals)
    #sigma_vals = np.array(sigma_vals)
    #timings = np.array(timings)
    #mem_uses = np.array(mem_uses)
    #energies = np.array(energies)
    madelungs = np.array(madelungs)

    #energies_kJmol = energies * N_A * 1e-3
    plt.plot(k_cut_vals, madelungs, 'x', label=r'$\sigma =$' + '{}'.format(np.around(sigma * 1e11, 2)) + r'$\cdot 10^{-11}$')

    if sigma == sigma_vals.max():
        plt.legend()
        plt.savefig(fname='../data/ewald_sum_sigma.png', dpi=900)
# %%