# Find most common fragments in a dataframe with SMILES column

The script finds the most common fragments from a given SMILES dataset that have more heavy atoms than a given threshold value. It first reads the dataset as a CSV file and extracts the SMILES column. Then it uses the BRICS algorithm to decompose each SMILES into fragments and removes any hydrogens present in the fragments. Finally, it counts the frequency of each fragment and outputs the n most common fragments along with their frequency. The output can also be saved as a CSV file and/or an SVG image showing the most common fragments. The script can be run from the command line with various arguments to control the output.

Depends on pandas and rdkit.

## Example usage

Find the most common 5 fragments which have more than 5 heavy atoms in the SMILES column in the test.csv:
```
common_frags.py -n 5 -m 5 -o 1 test.csv
```
