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
        self.particles, self.nuclei, \
            self.charges, self.coordinates, \
                self.N_part, self.volume, \
                    self.T_matrix, self.G_matrix \
                        = self._build_unit_cell()
            
    def _build_unit_cell(self):
        '''
        Initializes lists that contains particle names, nuclear charges, ionic charges and coordinates in the given crystal type.
        For now, only implemented for NaCl type crystals.
        '''
        if self.crystal_type == 'NaCl':
            particles, nuclei, charges, coordinates, N_part, volume, T_matrix, G_matrix = self._build_NaCl() 
        
        else: print('Cell types other than NaCl not yet implemented')

        return particles, nuclei, charges, coordinates, N_part, volume, T_matrix, G_matrix

    def _build_NaCl(self):
        '''
        Initializes lists that contains particle names, nuclear charges, ionic charges and coordinates
        for NaCl type crystals.

        See :ref [1]: Notes on Ewald summation techniques (docs).  
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

        # determine number of particles in the unit cell
        N_part = len(coordinates)

        # set up unit vectors to construct translation and reciprocal vectors
        x_vec = np.array([1, 0, 0])
        y_vec = np.array([0, 1, 0])
        z_vec = np.array([0, 0, 1])

        # set up translation vector matrix
        T_matrix = np.zeros((3, 3))
        
        # get individual components from ref [1]
        T_1 = 1/2 * (x_vec + y_vec)
        T_2 = 1/2 * (y_vec + z_vec)
        T_3 = 1/2 * (z_vec + x_vec)
        T_matrix[0, :] = T_1[:]
        T_matrix[1, :] = T_2[:]
        T_matrix[2, :] = T_3[:]

        # set up reciprocal vector matrix
        G_matrix = np.zeros((3, 3))

        # get individual components from ref [1]
        G_1 = 2 * np.pi * (x_vec + y_vec - z_vec)
        G_2 = 2 * np.pi * (-x_vec + y_vec + z_vec)
        G_3 = 2 * np.pi * (x_vec - y_vec + z_vec)
        G_matrix[0, :] = G_1[:]
        G_matrix[1, :] = G_2[:]
        G_matrix[2, :] = G_3[:]

        # determine unit cell volume
        volume = np.dot(T_1, np.cross(T_2, T_3))

        return particles, nuclei, charges, coordinates, N_part, volume, T_matrix, G_matrix

    def _build_cubic_lattice(self):
        '''
        Construct cubic periodic lattice from given cell parameter a. Each row is a lattice vector
        (a x a x a)
        '''
        lattice_matrix = np.eye(3, 3, dtype=np.float64) * self.cell_length

        inv_lat_matrix = np.linalg.inv(lattice_matrix)

        return lattice_matrix, inv_lat_matrix

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
        print(self.N_part)
        
        # initialize lists
        coords_in_cut = []
        charges_in_cut = []
        dist_in_cut = []

        # loop over all particles in all periodic images that may lie within cutoff radius
        for n_i in n:
            for j in range(self.N_part):
                r_j = self.coordinates[j] + n_i * self.cell_length
                r_ij = r_j - r_i
                d_ij = np.linalg.norm(r_ij)
                if d_ij > 0.0 and d_ij <= r_cut:
                    charges_in_cut.append(self.charges[j])
                    coords_in_cut.append(r_j)
                    dist_in_cut.append(d_ij)
        
        return charges_in_cut, coords_in_cut, dist_in_cut

    def get_T_G(self, n_cut, k_cut):
        '''
        Returns lists of translation and reciprocal repeat vector arrays for the given r_cut
        and k_cut values
        '''
        # get translation vector factors 
        n = self.get_n

        # initialize T_vectors list
        T_vectors = []

        # calculate T_vectors and append them to the list
        for n_i in n:
            T_vectors.append(n_i[0] * self.T_matrix[0] + n_i[1] * self.T_matrix[1] + n_i[2] * self.T_matrix[2])

        # get reciprocal repeat vector factors 
        k = self.get_k

        # initialize G_vectors list
        G_vectors = []

        # calculate G_vectors and append them to the list
        for k_i in k:
            G_vectors.append(k_i[0] * self.G_matrix[0] + k_i[1] * self.G_matrix[1] + k_i[2] * self.G_matrix[2])

        return T_vectors, G_vectors


    def _get_n(self, n_cut):
        '''
        Initializes list of arrays that contains all repeat vectors for a given n_cut.
        
        **Parameters:**
        n_cut : int
            Number that specifies the highest entry in any dimension that the repeat vectors may have.
            (Example: n_cut = 2 => n = [[-2, -2, -2], [-1, -2, -2], ..., [-1, -1, -1], ..., [0, 0, 0], [1, 0, 0], ..., [1, 1, 1]])
        '''
        n = []
        for i in np.arange(-n_cut, n_cut + 1):
            for j in np.arange(-n_cut, n_cut + 1):
                for k in np.arange(-n_cut, n_cut + 1):
                    n_ijk = np.array([i, j, k])   
                    n.append(n_ijk)        
        return n

    def _get_k(self, k_cut):
        '''
        Initializes list of arrays that contains all repeat vectors in reciprocal space for a given k_cut with k = [0, 0, 0] omitted.
        
        **Parameters:**
        k_cut : int
            Number that specifies the highest entry in any dimension that the repeat vectors may have.
            (Example: k_cut = 2 => n = [[-2, -2, -2], [-1, -2, -2], ..., [-1, -1, -1], ..., [1, 0, 0], ..., [1, 1, 1]])
        '''
        n = []
        for i in np.arange(-k_cut, k_cut + 1):
            for j in np.arange(-k_cut, k_cut + 1):
                for l in np.arange(-k_cut, k_cut + 1):
                    k_ijl = np.array([i, j, l])
                    
                    if np.linalg.norm(k_ijl) != 0:
                        n.append(k_ijl)        
        return n

# %%
