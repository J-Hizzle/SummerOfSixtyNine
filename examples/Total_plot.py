import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])
import time
import csv
import numpy as np
from scipy.constants import N_A
# %%
csvfile = '../data/direct_sum_results.csv'

# open the file in read mode
f = open(csvfile, 'r')

reader = csv.reader(csvfile)

# initialize parameter lists
n_cut_vals = []
timings = []
mem_uses = []
energies = []
madelungs = []

for line in reader:
    n_cut_vals.append(line[0])
    energies.append(line[1])
    madelungs.append(line[2])
    timings.append(line[3])
    mem_uses.append(line[4])

f.close()
# %%
n_cut_vals = np.array(n_cut_vals)
timings = np.array(timings)
mem_uses = np.array(mem_uses)
energies = np.array(energies)
madelungs = np.array(madelungs)

# %%
energies_kJmol = energies * N_A * 1e-3
print(energies_kJmol)
# %%
fig, ax = plt.subplots(2, 2, figsize=(10, 10))

plt.sca(ax[0, 0])
plt.plot(n_cut_vals[1:], energies_kJmol[1:], 'bo')
plt.ylabel(r'$E_{direct}$ in $\mathrm{kJ/mol}$')
plt.yticks((np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5)), np.around(np.linspace(energies_kJmol[1:].min(), energies_kJmol[1:].max(), 5), 2))
plt.title(r'Electrostatic energy')

plt.sca(ax[0, 1])
plt.plot(n_cut_vals[1:], madelungs[1:] , 'ro')
plt.yticks((np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5)), np.around(np.linspace(madelungs[1:].min(), madelungs[1:].max(), 5), 5))
plt.ylabel(r'$M_{direct}$')
plt.title(r'Madelung constant')

plt.sca(ax[1, 0])
plt.plot(n_cut_vals, timings, 'go')
plt.ylabel(r'$t_{run}$ in $\mathrm{s}$')
plt.title(r'Execution time')

plt.sca(ax[1, 1])
plt.plot(n_cut_vals, mem_uses, 'o', color='purple')
plt.ylabel(r'$m$ in $\mathrm{MiB}$')
plt.title(r'Memory usage')
plt.xlabel(r'$n_{cut}$')

fig.suptitle(r'$\mathbf{Direct sum results}$', y=0.94, fontsize=20)

fig.savefig(fname='../data/direct_sum_results', dpi=900)
# %%
