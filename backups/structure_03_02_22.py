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
        List of different species in the unit_cell.
    cell_length : float
        Characteristic length of the unit_cell.
    ions_per_cell : int
        Number of the ions of one species in the unit cell.
    particles : list, dtype = string, shape = (N_particles,)
        List of particle names in the unit_cell.
    nuclei : list, dtype = int, shape = (N_particles,)
        List of nuclear charges of the particles.
    charges : list, dtype = int, shape = (N_particles,)
        List of charges of the ions in the unit_cell.
    coordinates : list, dtype = tuple, shape = (N_particles,)
        List of ion coordinates in the unit_cell.
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

        # get coordinates from build_unit_cell function
        self.particles, self.nuclei, self.charges, self.coordinates = self.build_unit_cell()
            
    def _build_unit_cell(self):
        '''
        Initializes lists that contains particle names, nuclear charges, ionic charges and coordinates in the given crystal type.
        For now, only implemented for NaCl type crystals.
        '''
        coordinates = []

        if self.crystal_type == 'NaCl':
            # append list of non-reduntant coordinates for sodium ions
            coordinates.append([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])

            # append list of non-redundant coordinates for chloride ions
            coordinates.append([[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]])

        # set up lists with particle names, atomic numbers, charges
        particles = []
        nuclei = []
        charges = []
        
        for i in range(len(coordinates)):
            for j in range(len(coordinates[i])):
                particles.append(self.species[i])
                nuclei.append(self.atomic_numbers[i])
                charges.append(self.oxidation_states[i])

        # flatten coordinate list and scale according to unit_cell length
        coordinates = [np.asarray(coord_vec) * self.cell_length for spec_list in coordinates for coord_vec in spec_list]

        return particles, nuclei, charges, coordinates

    def _build_cubic_lattice(self):
        '''
        Construct cubic periodic lattice from given cell parameter a. Each row is a lattice vector
        (a x a x a)
        '''
        lattice_matrix = np.eye(3, 3, dtype=np.float64) * self.cell_length

        inv_lat_matrix = np.linalg.inv(lattice_matrix)

        return lattice_matrix, inv_lat_matrix


    def build_real_cut_spheres(self, real_cut, flat=True):
        '''
        Initializes list of arrays that contains all repeat vectors for a given cutoff value in real space.
        Obviously, one would like to limit the number of repeat vectors right from the start, but getting the correct set 
        is a bit tricky. At first, all repeat vectors for which the sum of the absolute values of the elements is smaller 
        or equal to the cutoff value are generated. Afterwards, the distances 
        For now, only implemented for n <= 1.
        '''
        cutoff_spheres = []

        if real_cut >= 0:
            cutoff_spheres.append([[0, 0, 0]])

        if real_cut >= 1:
            cutoff_spheres.append([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1], \
                                [-1, 0, 0], [0, -1, 0], [0, 0, -1], [-1, -1, 0], [-1, 0, -1], [0, -1, -1], [-1, -1, -1], \
                                [-1, 1, 0], [-1, 0, 1], [-1, 1, 1], \
                                [1, -1, 0], [0, -1, 1], [1, -1, 1], \
                                [1, 0, -1], [0, 1, -1], [1, 1, -1], \
                                [-1, -1, 1], [-1, 1, -1], [1, -1, -1]])
        
        # flatten the list and convert to arrays in proper units
        cutoff_spheres = [self.cell_length * np.asarray(rep_vec, dtype=np.float64) for cutoff_sphere in cutoff_spheres for rep_vec in cutoff_sphere]

        return cutoff_spheres

    def _get_n(n_cut):
        n = []
        for i in np.arange(-n_cut, n_cut):
            for j in np.arange(-n_cut, n_cut):
                for k in np.arange(-n_cut, n_cut):
                    n_ijk = np.array([i, j, k])
                    
                    if np.linalg.norm(n_ijk) != 0:
                        n.append(n_ijk)

        return n
# %%
