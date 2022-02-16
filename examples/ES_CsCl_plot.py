# %%
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import time
import csv
import pandas as pd
import numpy as np
from scipy.constants import N_A
# %%
csvfile = '../data/ewald_sum_CsCl_results.csv'

frame = pd.read_csv(csvfile, delimiter=',')

array = frame.to_numpy()

# get indices for respective sigma and n_cut values
n_cut = 2
sigma = 1e-10

n_cut_indices = array[:, 0] == n_cut
sigma_indices = array[:, 2] == sigma

indices = (array[:, 0] == n_cut) & (array[:, 2] == sigma)

# %%

#n_cut_vals = array[:, 0]
k_cut_vals = array[:, 1][indices]
#sigma_vals = array[:, 2]
energies = array[:, 3][indices]
madelungs = array[:, 4][indices]
timings = array[:, 5][indices]
mem_uses = array[:, 6][indices]
# %%
#n_cut_vals = np.array(n_cut_vals)
k_cut_vals = np.array(k_cut_vals)
#sigma_vals = np.array(sigma_vals)
timings = np.array(timings)
mem_uses = np.array(mem_uses)
energies = np.array(energies)
madelungs = np.array(madelungs)

energies_kJmol = energies * N_A * 1e-3
# %%
fig, ax = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True)

plt.sca(ax[0, 0])
plt.plot(k_cut_vals, energies_kJmol, 'bo')
plt.ylabel(r'$E_{ewald}$ in $\mathrm{kJ/mol}$')
plt.yticks((np.linspace(energies_kJmol.min(), energies_kJmol.max(), 5)), np.around(np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5), 2))
plt.title(r'Electrostatic energy')
plt.xlabel(r'$k_{cut}$')

plt.sca(ax[0, 1])
plt.plot(k_cut_vals, madelungs, 'ro')
plt.yticks((np.linspace(madelungs.min(), madelungs.max(), 5)), np.around(np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5), 5))
plt.ylabel(r'$M_{ewald}$')
plt.title(r'Madelung constant')
plt.xlabel(r'$k_{cut}$')

plt.sca(ax[1, 0])
plt.plot(k_cut_vals, timings, 'go')
plt.ylabel(r'$t_{run}$ in $\mathrm{s}$')
plt.title(r'Execution time')
plt.xlabel(r'$k_{cut}$')

plt.sca(ax[1, 1])
plt.plot(k_cut_vals, mem_uses, 'o', color='purple')
plt.ylabel(r'$m$ in $\mathrm{MiB}$')
plt.title(r'Memory usage')
plt.xlabel(r'$k_{cut}$')

fig.suptitle(r'Ewald sum results for CsCl ($\sigma = 5 \cdot 10^{-11}$, $n_{cut} = 2$)', fontsize=20, fontweight='bold')

fig.savefig(fname='../data/ewald_sum_results_CsCl_k_n{0}_s511.png'.format(n_cut, sigma), dpi=900)
# %%