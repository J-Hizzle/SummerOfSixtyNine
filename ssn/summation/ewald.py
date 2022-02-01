class EwaldSummation:
    '''
    class that gets data from structure class and calculates its Ewald sum,
    with the total electrostatic energy given as
    E_tot = E_real + E_recip + E_self + E_dip.

    For now, E_dip = 0.
    '''

    def __init__(self, structure, real_cutoff, recip_cutoff, alpha):
        '''
        

        
        '''