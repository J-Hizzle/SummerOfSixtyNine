'''
Main module of the SSN package that manages the user input and calls lower-level functions accordingly
'''
# %%
from SummerOfSixtyNine.ssn.summation.direct import direct_sum
from SummerOfSixtyNine.ssn.summation.ewald import ewald_sum
# %%

def run_ssn(struc, cutoff, n_cut, k_cut, sigma, technique, outval):
    '''
    Main function of the SSN package that calls lower-level functions.

    **Parameters:**
    
    outval : string, options=['E', 'F']
        Specifies return value of the function. Option include total energy (E) or force (F).
    technique : string, options=['DS', 'ES', 'PME']
        Specifies the technique to use for the calculation.
        Options include direct sum (DS), Ewald summation (ES), particle-mesh Ewald (PME).
    struc : Structure class instance
        Class containing all information about the struc of the given system.    
    '''
    # 1) input
    # first step make python script interface, maybe later with command line, but dunno
    # instantiate struc class

    # 2) calculation
    # implement these in lower-level modules within ssn directory
    
    # if technique='ES', instantiate Ewald summation class
    if technique == 'ES':
        E = ewald_sum(struc, sigma, n_cut, k_cut)

    if technique == 'DS':
        E = direct_sum(struc, cutoff)

    # 3) output
    # for now, put out results as return value
    if outval == 'E':
        return E