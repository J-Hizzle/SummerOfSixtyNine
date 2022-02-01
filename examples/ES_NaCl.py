"""
Example script to showcase the functionality of the implemented Ewald summation technique to calculate 
the electrostatic energy of the sodium cloride unit cell
"""
# %%
import SummerOfSixtyNine.ssn as ssn
# %%
# specify input parameters
outval      = 'E'       # calculate energy
technique   = 'ES'      # use Ewald summation

# instantiate structure class

# specify grid parameters


# initialize grid
# %%
# call main function with specified parameters
ssn.main(outval, technique, crystal)  
# %%
import csv
# %%
with open('../data/cell_parameters_dataset.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    for row in reader:
        print(row, '\n')
# %%
