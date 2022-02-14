'''
Main module of the SSN package that manages the user input and calls lower-level functions accordingly
'''
# %%
from SummerOfSixtyNine.ssn.summation.direct import direct_sum
from SummerOfSixtyNine.ssn.summation.ewald import ewald_sum
from SummerOfSixtyNine.ssn.madelung import madelung_constant
# %%
def run_ssn(structure, n_cut, technique, outval, k_cut=None, sigma=None):
    '''
    Main function of the SSN package that calls lower-level functions.

    **Parameters:**
    
    outval : list, dtype=string, options=['E', 'M']
        Specifies return value(s) of the function. Option include total energy (E) or Madelung constant (M).
    technique : string, options=['DS', 'ES', 'PME']
        Specifies the technique to use for the calculation.
        Options include direct sum (DS), Ewald summation (ES), particle-mesh Ewald (PME).
    structure : Structure class instance
        Class containing all information about the struc of the given system.    
    '''
    # 1) calculation    
    if technique == 'ES':
        # input error handling
        if sigma == None: print('Please provide Gaussian paramter for Ewald summation!')
        if k_cut == None: print('Please provide reciprocal cutoff for Ewald summation!')
        
        E = ewald_sum(structure, sigma, n_cut, k_cut)

    if technique == 'DS':
        E = direct_sum(structure, n_cut)

    # 2) output
    output = []

    if 'E' in outval:
        output.append(E)

    if 'M' in outval:
        M = madelung_constant(structure, E)
        output.append(M)

    return output