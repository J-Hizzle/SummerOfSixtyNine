'''
Translate fortran implementation from Ronald M. Pratt to python in order to understand
how setting up the crystal structure works.

See :ref [1]: R. M. Pratt, Madelung constants for ionic crystals using the Ewald sum, 2001.
'''
# %%
from math import *
import numpy as np
import scipy.constants as constants
# %%
# define calculation functions
def ewald_self(N_part, charges, alpha_d):
    E_self = 0.0

    for i in range(N_part):
        E_self += charges[i]**2

    E_self *= -alpha_d/np.sqrt(pi)
    return E_self

def ewald_real(N_part, coordinates, charges, r_cut, alpha_d):
    E_real = 0.0

    for i in range(N_part):
        for j in range(N_part):
            if i != j:
                # calculate distance r_ij
                r_ij = coordinates[j] - coordinates[i]

                # cutoff criterion
                r_ij_abs = np.linalg.norm(r_ij)

                if r_ij_abs < r_cut:
                    E_real += charges[i] * charges[j] * erfc(r_ij_abs * alpha_d)/r_ij_abs
    return E_real

def ewald_rec(N_part, kmax):
    # set up arrays and variables
    EXPIKR = np.zeros(N_part)
    E_l = np.zeros((N_part, 2 * kmax + 1), dtype=complex)
    E_m = np.zeros_like(E_l)
    E_n = np.zeros_like(E_l)
    TWOPI = 8 * atan(1.0)
    
    # set box size
    c_l, c_m, c_n = 2.0, 2.0, 2.0
    V = c_l * c_m * c_n

    E_rec = 0.0

    # store exponential factors
    for i in range(N_part):
        for k in np.arange(-kmax, kmax + 1):
            E_l[i, k] = complex(cos(k * TWOPI * coordinates[i][0]/c_l), sin(k * TWOPI * coordinates[i][0]))
            E_m[i, k] = complex(cos(k * TWOPI * coordinates[i][1]/c_l), sin(k * TWOPI * coordinates[i][1]))
            E_n[i, k] = complex(cos(k * TWOPI * coordinates[i][2]/c_l), sin(k * TWOPI * coordinates[i][2]))
    
    # loop over wave vectors
    for l in np.arange(-kmax, kmax + 1):
        r_l = TWOPI * l/c_l
        for m in np.arange(-kmax, kmax + 1):
            r_m = TWOPI * m/c_m
            for n in np.arange(-kmax, kmax + 1):
                r_n = TWOPI * n/c_n
                
                # test on magnitude of the k-vector
                kk = l**2 + m**2 + n**2 

                # skip if k vector is zero
                if kk != 0:
                    # calculate coefficient 
                    KSQ = r_l ** 2 + r_m ** 2 + r_n ** 2
                    AK = TWOPI/V * exp(-KSQ/(4 * alpha_d**2))/KSQ

                    # form exp(IKR) for each particle
                    for i in range(N_part):
                        EXPIKR[i] = E_l[i, l] * E_m[i, m] * E_n[i, n]

                    # form sums for each species
                    energy_sum = 0.0

                    for i in range(N_part):
                        energy_sum += charges[i] * EXPIKR[i]
                    # accumulate k-space potential energy
                    E_rec += AK * (energy_sum * energy_sum)

    return E_rec


# %%

Natoms  = 8
pi      = constants.pi
eps_0   = constants.epsilon_0
r_cut   = 0
alpha   = 0
facpe   = 0
kmax    = 26

# make array of non-redundant coordinates for ions
coordinates = np.asarray([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], \
                        [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]])

# make array of charges
charges = np.array([-1., -1., -1., -1., 1., 1., 1., 1.])

# number of atoms etc.
N_mol = 4
N_spec = 2
N_ion_1 = N_mol
N_ion_2 = N_mol
N_part = 2 * N_mol

# set dimensionless cutoff radius as half box length
r_cut = 0.5

# set charge distribution parameter
alpha = 8.0

# square of dimensionless cutoff
r_cut_square = r_cut**2

# dimensionless alpha 
alpha_d = alpha / 2

ewald_rec(N_part, kmax)
# %%

# calculate real-space contributions 
E_real = ewald_real()

# calculate rec-space contribution
E_rec = ewald_rec()

# calculate total energy 
E_tot = E_real + E_rec

# calculate madelung constant 
C_mad = E_tot/N_mol/charges[0]/charges[N_mol+1]












# %%
