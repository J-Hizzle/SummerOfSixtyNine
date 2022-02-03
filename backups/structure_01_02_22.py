# %%
import numpy as np
import re
# %%
class Structure:
    '''
    class that contains all information about crystal structure.

    **Attributes:**

    formula : string
        Chemical formula of the compound.
    atomic_numbers : list, dtype=int, shape = (N_species)
        Atomic numbers of each species given in the formula.
    oxidation_states : list, dtype=int, shape = (N_species)
        Charges of the respective species in the crystal.
    species : list, dtype = string, shape = (N_species)
        List of different species in the supercell.
    cell_length : float
        Characteristic length of the supercell.
    ions_per_cell : int
        Number of the ions of one species in the unit cell.
    particles : list, dtype = string, shape = (N_particles,)
        List of particle names in the supercell.
    nuclei : list, dtype = int, shape = (N_particles,)
        List of nuclear charges of the particles.
    charges : list, dtype = int, shape = (N_particles,)
        List of charges of the ions in the supercell.
    coordinates : list, dtype = tuple, shape = (N_particles,)
        List of ion coordinates in the supercell.
    '''

    def __init__(self, formula, atomic_numbers, oxidation_states, crystal_type, cell_parameters):
        '''
        Attach structure values to class. 

        **Parameters:**

        formula :           See above.
        oxidation_states :  See above.
        atomic_numbers :    See above.
        crystal_type : string, options = ['NaCl']
            Structure of the unit cell, necessary to specify the Bravais-lattice and coordinates.
        cell_parameters : list, dtype = (float, int), shape = (2,)
            First entry is the length of the unit cell L, second entry is coordination number Z.
        atom_positions : list, dtype = tuple, shape = (N_species)
            List of tuples where each tuple contains the position of the respective species in the unit cell.
        '''
        # initialize parameters directly from the input
        self.formula = formula
        self.oxidation_states = oxidation_states
        self.atomic_numbers = atomic_numbers
        self.cell_length = cell_parameters[0]
        self.ions_per_cell = cell_parameters[1]
        self.crystal_type = crystal_type

        # initialize species list from formula
        self.species = re.sub( r"([A-Z])", r" \1", self.formula).split()

        # get coordinates from build_supercell function
        self.particles, self.nuclei, self.charges, self.coordinates = build_supercell(self.species, self.atomic_numbers, self.oxidation_states, self.crystal_type, self.cell_length)
            
def build_supercell(species, atomic_numbers, oxidation_states, crystal_type, cell_length):
    '''
    Initializes lists that contains particle names, nuclear charges, ionic charges and coordinates in the given crystal type.
    For now, only implemented for NaCl type crystals.
    '''
    coordinates = []

    if crystal_type == 'NaCl':
        # append list of coordinates for sodium ions
        coordinates.append([[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0], [0, 1, 0], [1, 1, 0], \
                    [0.5, 0, 0.5], [0, 0.5, 0.5], [1, 0.5, 0.5], [0.5, 1, 0.5], \
                    [0, 0, 1], [1, 0, 1], [0.5, 0.5, 1], [0, 1, 1], [1, 1, 1]])

        # append list of coordinates for chloride ions
        coordinates.append([[0.5, 0, 0], [0, 0.5, 0], [1, 0.5, 0], [0.5, 1, 0], \
                            [0, 0, 0.5], [1, 0, 0.5], [0.5, 0.5, 0.5], [0, 1, 0.5], [1, 1, 0.5], \
                            [0.5, 0, 1], [0, 0.5, 1], [1, 0.5, 1], [0.5, 1, 1]])

    # set up lists with particle names, atomic numbers, charges
    particles = []
    nuclei = []
    charges = []
    
    for i in range(len(coordinates)):
        for j in range(len(coordinates[i])):
            particles.append(species[i])
            nuclei.append(atomic_numbers[i])
            charges.append(oxidation_states[i])

    # flatten coordinate list and scale according to supercell length
    coordinates = [np.asarray(coord_vec) * cell_length for spec_list in coordinates for coord_vec in spec_list]

    return particles, nuclei, charges, coordinates