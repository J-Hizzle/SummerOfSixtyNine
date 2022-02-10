from scipy.constants import epsilon_0, pi, e
from scipy.special import erfc
import numpy as np
from math import exp, cos, sin

def ewald_sum(structure, sigma, real_cut, rec_cut):
    '''
    Calculates electrostatic energy for a given structure using the ewald summation technique
    according to ref [1]

    **Parameters:**

    structure : Structure class instance
        Object containing all information about the structure of the given system.
    sigma : float
        Gaussian standard deviation. 
    real_cut : float
        Cutoff radius in real space.
    rec_cut : float
        Cutoff for distance in reciprocal space.
    
    See :ref [1]: Notes on Ewald summation techniques (docs).
    '''
    # set summation parameters
    q = structure.charges
    r = structure.coordinates
    
    # set number of particles in unit cell
    N = len(q)

    # get repeat vectors in real space T
    T_vectors, G_vectors = structure.get_T_G(real_cut, rec_cut)

    # initialize total electrostatic energy 
    E_tot = 0.0

    # loop over particles in unit cell
    for i in range(N):
        for j in range(N):
            # calculate self-energy term 
            E_self = ewald_self(sigma, i, j)

            # calculate short-range term in real space
            E_short = ewald_short(structure, sigma, T_vectors, i, j)
            
            # calculate long-range term in reciprocal space
            E_long = ewald_long(structure, sigma, G_vectors, i, j)
            
            # add sum over the contributions to total energy
            E_tot += E_long + E_short - E_self

    # calculate units and compensate for double-counting
    E_tot *= e**2/(4 * pi * epsilon_0) * 1/2

    return E_tot

def ewald_self(sigma, i, j):
    '''
    Returns self-energy of i-th and j-th particle
    '''
    if i == j:
        E_self = 1/np.sqrt(sigma * pi)
    else:
        E_self = 0.0
    
    return E_self

def ewald_short(structure, sigma, T_vectors, i, j):
    '''
    Calculates real space contribution to electrostatic energy, consisting of the short-range part 
    and the self-energy term.
    '''
    # set summation parameters
    r = structure.coordinates
    
    # initialize real space energy
    E_short = 0.0
    
    for T_i in T_vectors:
        # calculate only if it's not the same particle
        if not (np.linalg.norm(T_i) == 0 and i == j):
            d_ij = np.linalg.norm(r[j] - r[i] + T_i)
            E_short += erfc(1/(2 * np.sqrt(sigma)) * d_ij)/d_ij
            
    return E_short

def ewald_long(structure, sigma, G_vectors, i, j):
    '''
    Calculates reciprocal space contribution to the electrostatic energy.
    '''
    # set summation parameters
    r = structure.coordinates
    V = structure.volume
    
    # initialize reciprocal space energy
    E_long = 0.0

    for G_i in G_vectors:
        # omit k = 0
        if not (np.linalg.norm(G_i) == 0):
            r_ij = r[j] - r[i]
            k_abs = np.linalg.norm(G_i)
            E_long += complex(cos(G_i * r_ij), sin(G_i * r_ij)) * exp(-k_abs**2 * sigma)/k_abs**2

    # include prefactor
    E_long *= 4 * pi/V

    return E_long