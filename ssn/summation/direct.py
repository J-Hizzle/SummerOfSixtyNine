# %%
from scipy.constants import epsilon_0, pi, e
import numpy as np
# %%
def direct_sum(structure, n_cut):
    '''
    Calculates electrostatic energy for a given structure using the direct sum method
    with real_cut according to

    **Parameters:**

    structure : Structure class instance
        Class containing all information about the structure of the given system.
    n_cut : float
        Cubic cutoff shell in real space.
    
    See :ref [1]: H. B. Lee, W. Cai, Ewald Summation for Coulomb Interactions in a Periodic Supercell, 2009.
    '''
    # set summation parameters
    L = structure.cell_length
    q = structure.charges
    r = structure.coordinates
    N = structure.N_part
    
    # set repeat vectors
    n = structure.get_rep(n_cut)

    # initialize electrostatic energy
    E = 0.0
    
    # add direct summands succeedingly
    for i in range(N):
        for j in range(N):
            for n_i in n:
                nL = n_i * L
                d_ij = np.linalg.norm(r[j] - r[i] + nL)
                if not (d_ij == 0.0):
                    E += q[i] * q[j] / (d_ij)
    
    # account for double counting and multiply with constants
    E *= e**2/(4 * pi * epsilon_0) * 1/2

    return E