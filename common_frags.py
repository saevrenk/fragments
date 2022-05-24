"""
Find the most common fragments for a given SMILES set

Given a dataset in .csv format with SMILES column, returns 
the n most common fragments which have more heavy atoms than m 

"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import BRICS
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)  # set to avaoid print of warnings


def smiles2commonfrags(smiles, m):
    """
    get fragments with more than m heavy atoms and their frequency for a smiles list

    args:
       smiles - a list of SMILES
       m      - heavy atom numbers in a fragment

    returns:
       dictionary of fragments larger than m and their frequency

    """
    lfrag_collect = {}

    dummy = Chem.MolFromSmiles("*")

    for smile in smiles:

        mol = Chem.MolFromSmiles(smile)
        frags = [x for x in BRICS.BRICSDecompose(mol)]

        for f in frags:

            frag = AllChem.ReplaceSubstructs(
                Chem.MolFromSmiles(f), dummy, Chem.MolFromSmiles("[H]"), True
            )[0]
            frag = Chem.RemoveHs(frag)
            frag_smi = Chem.MolToSmiles(frag)

            if frag.GetNumHeavyAtoms() > m:
                if frag_smi in lfrag_collect.keys():
                    lfrag_collect[frag_smi] = lfrag_collect[frag_smi] + 1
                else:
                    lfrag_collect[frag_smi] = 1

    return lfrag_collect




if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(usage="%(prog)s [options] dataset.csv")
    parser.add_argument("dataset", metavar="dataset", type=str)
    parser.add_argument("-n", "--nfrags", help="n most common fragments")
    parser.add_argument("-m", "--minheavyatoms", help="smallest frag size to return")
    parser.add_argument(
        "-o", "--output_png", help="1/0 to output svg file showing fragments"
    )

    args = parser.parse_args()

    # read dataset
    df = pd.read_csv(args.dataset)
  
    if args.minheavyatoms is not None:
        m = int(args.minheavyatoms)
    else:
        m = 5

    if args.nfrags is not None:
        n = int(args.nfrags)
    else:
        n = 10

    if args.output_png is not None:
        o = int(args.output_png)
    else:
        o = 0
 
    # get dictionary of fragments
    frag_dict = smiles2commonfrags(df.SMILES, m)

    # sort dictionary to get the most common fragments
    sorted_frags = sorted(frag_dict.items(), key=lambda x: x[1], reverse=True)
    highfreq_smi = []
    freqs = []
    for (frag, freq) in sorted_frags[0:n]:
        freqs.append(freq)
        highfreq_smi.append(frag)

    # Convert fragment dictionary to dataframe and save as .csv file:
    df_frags = pd.DataFrame(
        list(frag_dict.items()), columns=["Frag_Smiles", "Frag_Frequency"]
    )
    df_frags.to_csv("fragment_freq.csv")

    if o == 1:

        frag_m = [Chem.MolFromSmiles(x) for x in highfreq_smi]
        svg = Draw.MolsToGridImage(
            frag_m, molsPerRow=4, useSVG=True, legends=[str(x) for x in freqs]
        )

        print("here")
        with open("frags.svg", "w") as figure:
            figure.write(svg)

    if len(frag) < n:
        print("\nOnly found %d fragments larger than %d heavy atoms!:\n" % (len(frag),m))
    else:
        print("\n%d most common fragments that have more heavy atoms than %d:\n" % (n, m))

    print("Fragment SMILE, frequency")
    [print("%s: %d" % (fr, fq)) for fr, fq in sorted_frags[0:n]]
