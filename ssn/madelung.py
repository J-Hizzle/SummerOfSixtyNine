from scipy.constants import epsilon_0, pi, e, N_A

def madelung_constant(structure, E):
    '''
    Calculate Madelung constant for a given structure and energy to validate the result.

    structure : Structure class instance
        Class containing all information about the structure of the given system.
    E : float
        Total electrostatic energy calculated for the given structure.
    '''
    # set parameters for formulas
    L = structure.cell_length
    d = structure.min_dist
    Z = structure.ions_per_cell
    k = 1/(4 * pi * epsilon_0)

    # calculate Madelung constant
    M = E * L * d/(e**2 * k * Z)

    return M