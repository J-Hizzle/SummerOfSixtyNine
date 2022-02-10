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
    e = constants.e

    # set summation parameters
    L = structure.cell_length
    q = structure.charges
    r = structure.coordinates
    
    # set repeat vectors
    n_cut = int(np.rint(real_cut/L + 1))
    n = structure.get_n(n_cut)

    # set number of particles in unit cell
    N = len(q)

    # initialize electrostatic energy
    E = 0.0
    
    # add direct summands succeedingly
    for i in range(N):
        charges_in_cut, coords_in_cut, dist_in_cut = structure.get_particles_in_cut(r[i], real_cut)
        for j in len(coords_in_cut):
            E += q[i] * charges_in_cut[j] / (dist_in_cut[j])
    
    # account for double counting and multiply with constants
    E *= e**2/(4 * pi * eps_0) * 1/2

    return E
# %%
