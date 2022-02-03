# %%
import scipy.constants as constants
import numpy as np
# %%
def direct_sum(structure, real_cut):
    '''
    Calculates electrostatic energy for a given structure using the direct sum method
    with real_cut according to

    **Parameters:**

    structure : Structure class instance
        Class containing all information about the structure of the given system.
    real_cut : float
        Cutoff radius in real space.
    
    See :ref [1]: H. B. Lee, W. Cai, Ewald Summation for Coulomb Interactions in a Periodic Supercell, 2009.
    '''
    # set physical constants
    eps_0 = constants.epsilon_0
    pi = constants.pi

    # set summation parameters
    L = structure.cell_length
    q = structure.charges
    r = structure.coordinates

    # set repeat vectors
    n = structure.build_real_cut_spheres(real_cut)

    # set number of particles within real_cut
    N = len(q) * len(n)

    # initialize electrostatic energy
    E = 0
    
    # add direct summands succeedingly
    for n_i in n:
        for i in range(N):
            for j in range(N):
                if i != j:
                    r_ij = r[j] - r[i]
                    nL = n_i * L
                    #print('r_ij({0},{1}) ='.format(i, j), r_ij)
                    #print('nL({0},{1}) ='.format(i, j), nL)
                    E += q[i] * q[j] / (np.linalg.norm(r_ij + nL))

    # account for double counting and multiply with constants
    E *= 1/(4 * pi * eps_0) * 1/2

    return E
# %%
