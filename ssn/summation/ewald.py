from scipy.constants import epsilon_0, pi, e
from scipy.special import erfc
import numpy as np

def ewald_sum(structure, sigma, n_cut, k_cut):
    '''
    Calculates electrostatic energy for a given structure using the ewald summation technique
    according to ref [1]

    **Parameters:**

    structure : Structure class instance
        Object containing all information about the structure of the given system.
    sigma : float
        Gaussian standard deviation. 
    n_cut : float
        Cubic cutoff shell in real space.
    k_cut : float
        Cubic cutoff shell in reciprocal space.
    
    See :ref [1]: Notes on Ewald summation techniques (docs).
        :ref [2]: Report on Ewald summation (docs).
    '''
    # set summation parameters
    r = structure.coordinates
    q = structure.charges
    N = structure.N_part
    L = structure.cell_length
    V = structure.volume

    # get real and reciprocal space repeat vectors
    n = structure.get_rep(n_cut)
    k = structure.get_rep(k_cut)

    # initialize total energy
    E_tot = 0.0

    # debugging
    E_long_tot = 0.0

    # loop over particles
    for i in range(N):
        for j in range(N):
            # calculate different contributions to energy
            E_long = ewald_long(r[i], r[j], k, L, V, sigma)
            E_short = ewald_short(r[i], r[j], n, L, sigma)
            E_self = ewald_self(i, j, sigma)

            # sum over all individual contributions and multiply with product of charges
            E_tot += q[i] * q[j] * (E_long + E_short - E_self)
            E_long_tot += E_long

    # calculate units and compensate for double-counting
    E_tot *= e**2/(8 * pi * epsilon_0)

    print('E_long_tot =', E_long_tot)

    return E_tot

def ewald_long(r_i, r_j, k, L, V, sigma):
    '''
    Calculates long-range interaction energy in reciprocal space.
    '''
    # initialize E_long
    E_long = 0.0 

    # loop over k-vectors
    for k_i in k:
        # calculate correct dimension of k
        k_dim = 2 * pi/L * k_i

        # calculate absolute value of k
        k_abs = np.linalg.norm(k_dim)

        # omitt k = 0 vector 
        if not (k_abs == 0.0):
            # add energy term to long-range energy
            E_long += np.cos(np.dot(k_dim, (r_j - r_i))) * np.exp(-sigma**2 * k_abs**2/2)/k_abs**2
            #print('E_long({0}) ='.format(k_i), E_long)
            #print('cos_term({0}) ='.format(k_i), np.cos(np.dot(k_dim, (r_j - r_i))))
            #print('factor_term({0}) ='.format(k_i), np.exp(-sigma**2 * k_abs**2/2)/k_abs**2)
            #print('exponent({0}) ='.format(k_i), -sigma**2 * k_abs**2/2)

    # account for prefactor
    E_long *= 4 * pi/V

    #print('E_long =', E_long)

    return E_long

def ewald_short(r_i, r_j, n, L, sigma):
    '''
    Calculates short-range interaction energy in real space.
    '''
    # initialize E_short
    E_short = 0.0

    # loop over n-vectors
    for n_i in n:
        # calculate correct dimensions of n
        nL = n_i * L

        # calculate absolute distance between particles
        d_ij = np.linalg.norm(r_j - r_i + nL)

        # omitt self-interaction
        if not (d_ij == 0.0):
            E_short += erfc(d_ij/(np.sqrt(2) * sigma))/d_ij
    
    #print('E_short =', E_short)
    
    return E_short

def ewald_self(i, j, sigma):
    '''
    Calculates self-energy term.
    '''
    if i == j:
        E_self = np.sqrt(2/pi)/sigma
    else:
        E_self = 0.0

    #print('E_self =', E_self)

    return E_self