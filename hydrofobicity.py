import pandas as pd
from biopandas.pdb import PandasPdb
from Bio import SeqIO
import os


def hydro_fob(file, type, scale):
    if scale=="Kyte-Doolittle":
        scale = {
            'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
            'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
            'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
            'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
        sekw = list(map(lambda x: list(str(x.seq)), SeqIO.parse(file, type)))
        furn = list(map(lambda x: (list(map(lambda n: scale[n], x))), sekw))

        return furn

    else:
        return AttributeError

def hydrofobicity_pdb(file):
    pdb = PandasPdb().read_pdb(file)
    df = pdb.df["ATOM"]
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    print(df.head())
    #residue_name for names of aminoacids and b_factor will be changed to hydrofobicity in relation to AM
    scale_pdb = {
        'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
        'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
        'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
        'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
    }
    sekw = list(map(lambda x: scale_pdb[x], df["residue_name"]))
    df["b_factor"] = sekw
    pdb.df["ATOM"] = df
    return pdb




if __name__ == '__main__':
    place = "/home/arjissuan/Desktop/tavle_significant/pdb"
    #print((hydrofobicity_pdb("/home/arjissuan/Desktop/tavle_significant/pdb/7715pdb/7715.pdb")))
    for folder in os.scandir(place):
        if folder.is_dir():
            for file in os.listdir(os.path.join(place, folder)):
                hydrofobicity_pdb(os.path.join(place, folder,file)).to_pdb(os.path.join(place,folder,"{}_new_factor.pdb".format(file.split(".")[0])))





    #anno_names = list("{}".format(item.split(".")[0]) for item in os.listdir("/home/arjissuan/Desktop/tavle_significant/Matrixes_with_statistics"))
    #for item in anno_names:
    #   print(len(hydro_fob("/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/{}_a_e_l_s_bez.fasta".format(item), 'fasta', "Kyte-Doolittle")))



