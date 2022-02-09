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
        self.particles, self.nuclei, self.charges, self.coordinates = self._build_unit_cell()
            
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


    def get_n(self, n_cut):
        '''
        Initializes list of arrays that contains all repeat vectors for a given n_cut with n = [0, 0, 0] omitted.
        
        **Parameters:**
        n_cut : int
            Number that specifies the highest entry in any dimension that the repeat vectors may have.
            (Example: n_cut = 2 => n = [[-2, -2, -2], [-1, -2, -2], ..., [-1, -1, -1], ..., [1, 0, 0], ..., [1, 1, 1]])
        '''
        n = []
        for i in np.arange(-n_cut, n_cut):
            for j in np.arange(-n_cut, n_cut):
                for k in np.arange(-n_cut, n_cut):
                    n_ijk = np.array([i, j, k])
                    
                    if np.linalg.norm(n_ijk) != 0:
                        n.append(n_ijk)        
        return n

    def get_particles_in_cut(self, r_i, r_cut):
        '''
        Returns list of charges, list of coordinate-arrays and list of distances of particles that lie within a distance r_cut to 
        a given particle at position r_i. 

        **Parameters:**

        r_i : np.ndarray, dtype = float64, shape = (3,)
            Position of a given particle.
        r_cut : float
            Cutoff radius to which particles will be considered.
        '''
        # set repeat vectors
        n_cut = int(np.rint(r_cut/self.cell_length + 1))
        n = self.get_n(n_cut)
        
        # initialize lists
        coords_in_cut = []
        charges_in_cut = []
        dist_in_cut = []

        # loop over all particles in all periodic images that may lie within cutoff radius
        for n_i in n:
            for j in len(self.coordinates):
                r_j = self.coordinates[j] + n_i
                r_ij = r_j - r_i
                d_ij = np.linalg.norm(r_ij)
                if d_ij > 0.0 and d_ij <= r_cut:
                    charges_in_cut.append(self.charges[j])
                    coords_in_cut.append(r_j)
                    dist_in_cut.append(d_ij)
        
        return charges_in_cut, coords_in_cut, dist_in_cut
# %%
