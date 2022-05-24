# Find most common fragments in a dataframe with SMILES column

Given a dataset in the form of a .csv file, the code returns n most common fragments having more than  m heavy atoms.  

Depends on pandas and rdkit.

## Example usage

Find the most common 5 fragments which have more than 5 heavy atoms in the dataframe test.csv SMILES column:

common_frags.py -n 5 -m 5 -o 1 test.csv
 
