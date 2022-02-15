import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import time
import csv
import pandas as pd
import numpy as np
from scipy.constants import N_A
# %%
csvfile = '../data/direct_sum_results.csv'

frame = pd.read_csv(csvfile, delimiter=',')

array = frame.to_numpy()

n_cut_vals = array[:, 0]
energies = array[:, 1]
madelungs = array[:, 2]
timings = array[:, 3]
mem_uses = array[:, 4]

# %%
n_cut_vals = np.array(n_cut_vals)
timings = np.array(timings)
mem_uses = np.array(mem_uses)
energies = np.array(energies)
madelungs = np.array(madelungs)

energies_kJmol = energies * N_A * 1e-3
# %%
fig, ax = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True)

plt.sca(ax[0, 0])
plt.plot(n_cut_vals[1:], energies_kJmol[1:], 'bo')
plt.ylabel(r'$E_{direct}$ in $\mathrm{kJ/mol}$')
plt.yticks((np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5)), np.around(np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5), 2))
plt.title(r'Electrostatic energy')
plt.xlabel(r'$n_{cut}$')

plt.sca(ax[0, 1])
plt.plot(n_cut_vals[1:], madelungs[1:] , 'ro')
plt.yticks((np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5)), np.around(np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5), 5))
plt.ylabel(r'$M_{direct}$')
plt.title(r'Madelung constant')
plt.xlabel(r'$n_{cut}$')

plt.sca(ax[1, 0])
plt.plot(n_cut_vals, timings, 'go')
plt.ylabel(r'$t_{run}$ in $\mathrm{s}$')
plt.title(r'Execution time')
plt.xlabel(r'$n_{cut}$')

plt.sca(ax[1, 1])
plt.plot(n_cut_vals, mem_uses, 'o', color='purple')
plt.ylabel(r'$m$ in $\mathrm{MiB}$')
plt.title(r'Memory usage')
plt.xlabel(r'$n_{cut}$')

fig.suptitle(r'Direct sum results', fontsize=20, fontweight='bold')

fig.savefig(fname='../data/direct_sum_results', dpi=900)
# %%