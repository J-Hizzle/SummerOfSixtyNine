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
    crystal_type : string
        Type of the crystal system
    N_part : int
        Number of particles in the unit cell.
    particles : list, dtype = string, shape = (N_particles,)
        List of particle names in the unit cell.
    volume : float
        Unit cell volume.
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
        atomic_numbers :    See above.
        oxidation_states :  See above.
        crystal_type : string, options = ['NaCl']
            Structure of the unit cell, necessary to specify the Bravais-lattice and coordinates.
        cell_parameters : list, dtype = (float, int), shape = (2,)
            First entry is the length of the unit cell L, second entry is coordination number Z.        '''
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
                self.N_part, self.volume \
                        = self._build_unit_cell()
            
    def _build_unit_cell(self):
        '''
        Initializes lists that contains particle names, nuclear charges, ionic charges and coordinates in the given crystal type.
        For now, only implemented for NaCl type crystals.
        '''
        if self.crystal_type == 'NaCl':
            particles, nuclei, charges, coordinates, N_part, volume = self._build_NaCl() 
        
        else: print('Cell types other than NaCl not yet implemented')

        return particles, nuclei, charges, coordinates, N_part, volume

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

        # determine unit cell volume
        volume = self.cell_length**3

        return particles, nuclei, charges, coordinates, N_part, volume

    def get_rep(self, cutoff):
        '''
        Initializes list of arrays that contains all repeat vectors for a given cutoff value.
        
        **Parameters:**
        cutoff : int
            Number that specifies the highest entry in any dimension that the repeat vectors may have.
            (Example: n_cut = 2 => n = [[-2, -2, -2], [-1, -2, -2], ..., [-1, -1, -1], ..., [0, 0, 0], [1, 0, 0], ..., [1, 1, 1]])
        '''
        n = []
        for i in np.arange(-cutoff, cutoff + 1):
            for j in np.arange(-cutoff, cutoff + 1):
                for k in np.arange(-cutoff, cutoff + 1):
                    n_ijk = np.array([i, j, k])   
                    n.append(n_ijk)        
        return n
# %%
