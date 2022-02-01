# %%
import scipy.constants as constants
import numpy as np
# %%
def direct_sum(structure, cutoff):
    '''
    Calculates electrostatic energy for a given structure using the direct sum method
    with cutoff according to

    **Parameters:**

    structure : Structure class instance
        Class containing all information about the structure of the given system.
    cutoff : float
        Number of spheres to consider in the direct sum around the unit cell.
    
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
    n = build_cutoff_spheres(cutoff)

    # set number of particles within cutoff
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

def build_cutoff_spheres(cutoff, flat=True):
    '''
    Initializes nested list of tuples that contains all repeat vectors in each cutoff-sphere for a given cutoff value
    and flattens it if requested.
    For now, only implemented for n <= 1.
    '''
    cutoff_spheres = []

    if cutoff >= 0:
        cutoff_spheres.append([[0, 0, 0]])

    if cutoff >= 1:
        cutoff_spheres.append([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1], \
                            [-1, 0, 0], [0, -1, 0], [0, 0, -1], [-1, -1, 0], [-1, 0, -1], [0, -1, -1], [-1, -1, -1], \
                            [-1, 1, 0], [-1, 0, 1], [-1, 1, 1], \
                            [1, -1, 0], [0, -1, 1], [1, -1, 1], \
                            [1, 0, -1], [0, 1, -1], [1, 1, -1], \
                            [-1, -1, 1], [-1, 1, -1], [1, -1, -1]])
    
    # flatten the list if requested
    if flat == True:
        cutoff_spheres = [np.asarray(rep_vec) for cutoff_sphere in cutoff_spheres for rep_vec in cutoff_sphere]

    return cutoff_spheres
# %%
