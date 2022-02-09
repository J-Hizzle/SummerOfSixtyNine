from scipy.constants import epsilon_0, pi, e
from scipy.special import erf, erfc
import numpy as np

from SummerOfSixtyNine.ssn.pratt import E_real

def ewald(structure, sigma, real_cut, rec_cut):
    '''
    Calculates electrostatic energy for a given structure using the direct sum method
    with real_cut according to

    **Parameters:**

    structure : Structure class instance
        Object containing all information about the structure of the given system.
    sigma : float
        Gaussian standard deviation. 
    real_cut : float
        Cutoff radius in real space.
    rec_cut : float
        Cutoff for distance in reciprocal space.
    
    See :ref [1]: H. B. Lee, W. Cai, Ewald Summation for Coulomb Interactions in a Periodic Supercell, 2009.
    '''
    # calculate real space contribution to electrostatic energy
    E_real = ewald_real(structure, sigma, real_cut)

    # calculate reciprocal space contribution to electrostatic energy 
    E_rec = ewald_rec(structure, sigma, rec_cut)

    # calculate total electrostatic energy
    E = E_real + E_rec

    return E

def ewald_real(structure, sigma, real_cut):
    '''
    Calculates real space contribution to electrostatic energy, consisting of the short-range part 
    and the self-energy term.
    '''
    # set summation parameters
    L = structure.cell_length
    q = structure.charges
    r = structure.coordinates
    
    # set repeat vectors
    n_cut = int(np.rint(real_cut/L + 1))
    n = structure.get_n(n_cut)

    # set number of particles in unit cell
    N = len(q)

    # initialize real space energy
    E_real = 0.0
    
    # loop over particles and calculate particles within cutoff sphere on the go
    for n_i in n:
        for i in range(N):
            for j in range(N):
                r_ij = r[j] - r[i]
                d = np.linalg.norm(r_ij + nL)
                # calculate self-energy
                if d == 0.0:
                    E_real += 1/(np.sqrt(2 * pi) * sigma) * q[i]**2

                # calculate short-range terms
                if i != j:
                    if d <= real_cut:
                        nL = n_i * L
                        
                        E_real += q[i] * q[j] / d * erfc(d/(np.sqrt(2) * sigma))

    # account for double counting and multiply with constants
    E_real *= e**2/(4 * pi * epsilon_0) * 1/2

    return E_real

def ewald_rec(structure, sigma, rec_cut):
    '''
    Calculates reciprocal space contribution to the electrostatic energy.
    '''
