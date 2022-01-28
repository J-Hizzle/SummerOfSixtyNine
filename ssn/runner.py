'''
Main module of the SSN package that manages the user input and calls lower-level functions accordingly
'''

def main(outval, technique, crystal):
    '''
    Main function of the SSN package that calls lower-level functions.

    **Parameters:**
    
    outval : string, options=['E', 'F']
        Specifies return value of the function. Option include total energy (E) or force (F).
    technique : string, options=['DS', 'ES', 'PME']
        Specifies the technique to use for the calculation.
        Options include direct sum (DS), Ewald summation (ES), particle-mesh Ewald (PME).
    crystal : string, options=['NaCl']
        Specifies the crystal for which to calculate the given output value with the given technique.
        Options inclue sodium chloride (NaCl).
    '''
    # 1) input
    # first step make python script interface, maybe later with command line, but dunno

    # 2) calculation
    # implement these in lower-level modules within ssn directory

    # 3) output
    # for now, put out results as return value
    return