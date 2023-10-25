from Bio import SeqIO, SearchIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import mdtraj as md

class MarrBacines:
    def __init__(self, location="/home/arjissuan/Desktop/", probes=("DBL", "IK1", "M11", "MD19A", "MD87", "MD89A", "MKC", "MSI", "SD1", "D7", "AAT")):
        self.probes = probes
        self.location = location
        self.Probe_name = "PROKKA_all.faa"
        self.DB_name = "mono_out.fsa"

    def DB_prep_for_phmmer(self):
        #probably main
        pass

    def find_files(self, probe, results_name): #now will give only one file
        if "*" in results_name:
            results_name = results_name.strip("*")
        files = os.listdir(os.path.join(self.location, probe))
        name_of_res = lambda x: x if results_name in x else None
        name = list(map(name_of_res, files))
        while None in name:
            name.remove(None)
        return os.path.join(self.location, probe, name[0])

    def domtab_to_excel(self, *results_name):
        #kon_choices
        funk = lambda x: (x.name, x.seq)
        for probe in self.probes:
            if results_name == ():
                file = os.path.join(self.location, probe, f"results_{probe}.txt")
                AMP_DB = dict(map(funk, SeqIO.parse(os.path.join(self.location, probe, self.DB_name), "fasta")))
            else:
                file = self.find_files(probe, "results_"+results_name[0])
                AMP_DB = dict(map(funk, SeqIO.parse(self.find_files(probe, results_name[1]), "fasta")))

            results = list(map(lambda x: x, SearchIO.parse(file, 'phmmer3-domtab')))
            df = pd.DataFrame(columns=["query_id", "bit-score", "bias", "e-value", "hit_id", "hit_description"])
            query_id = []
            bit = []
            bias = []
            evalue = []
            hit_id = []
            hit_descr = []
            for i in range(len(results)):
                for j in range(len(results[i])):
                    query_id.append(results[i][j].query_id)
                    bit.append(results[i][j].bitscore)
                    bias.append(results[i][j].bias)
                    evalue.append(results[i][j].evalue)
                    hit_id.append(results[i][j].id)
                    hit_descr.append(results[i][j].description)
            df["query_id"] = query_id
            df["bit-score"] = bit
            df["bias"] = bias
            df["e-value"] = evalue
            df["hit_id"] = hit_id
            df["hit_description"] = hit_descr
            AMP_probe = dict(map(funk, SeqIO.parse(os.path.join(self.location, probe, self.Probe_name), "fasta")))
            AMP_DB_seq = []
            AMP_probe_seq = []
            for h_id in df["hit_id"]:
                if h_id in AMP_DB.keys():
                    AMP_DB_seq.append(AMP_DB[h_id])
            for q_id in df["query_id"]:
                if q_id in AMP_probe.keys():
                    AMP_probe_seq.append(AMP_probe[q_id])
            df["query_seq"] = AMP_probe_seq
            df["hit_seq"] = AMP_DB_seq

            if results_name == ():
                df.to_excel(os.path.join(self.location, probe, f"out_{probe}.ods"))

            else:
                df.to_excel(os.path.join(self.location, probe, f"out_{probe}_{results_name[0]}.ods"))


    def AMP_DB_chem_prop(self,DB_name):
        AMP_DB = list(map(lambda x: x, SeqIO.parse(self.find_files(self.probes[0], DB_name), "fasta")))
        sec_lam = lambda x: {"sequence": (str(x.seq)),
                             "id": (x.name),
                             'length': (len(x.seq)),
                             "aromaticity": (ProteinAnalysis(str(x.seq)).aromaticity()),
                             "hydrofobicity": (ProteinAnalysis(str(x.seq)).gravy()),
                             "instability_index": (ProteinAnalysis(str(x.seq)).instability_index()),
                             # stable <40 unstable >40
                             "masa": (ProteinAnalysis(str(x.seq)).molecular_weight()),
                             "sec_structure_frac": (ProteinAnalysis(str(x.seq)).secondary_structure_fraction()),
                             "flex": (ProteinAnalysis(str(x.seq)).flexibility()),
                             "izoel_point": (IsoelectricPoint(x.seq).pi()),
                             }
        AMP_properties = list(map(sec_lam, AMP_DB))
        return pd.DataFrame(AMP_properties)

    def AMP_chem_prop(self, probe, name):
        value_prop = lambda x, y: {"sequence": x,
                                   "AMP_id": y,
                                   "length": len(x),
                                   'aromaticity': ProteinAnalysis(x).aromaticity(),
                                   'hydrofobicity': ProteinAnalysis(x).gravy(),
                                   'instability_index': ProteinAnalysis(x).instability_index(),
                                   'mass': ProteinAnalysis(x).molecular_weight(),
                                   'sec_structure_frac': ProteinAnalysis(x).secondary_structure_fraction(),
                                   'flexibility': ProteinAnalysis(x).flexibility(),
                                   'izoel_point': IsoelectricPoint(x).pi()}
        Df = pd.read_excel(self.find_files(probe, name), engine="odf")
        new_df = list(map(value_prop, Df["query_seq"], Df["query_id"]))
        return pd.DataFrame(new_df)

    def absolute(self, probe, out_name):
        Df = pd.read_excel(self.find_files(probe, out_name), engine="odf")
        print(probe, out_name)
        query_seq = np.array(list(map(lambda x: [len(x)], Df["query_seq"])))
        hit_seq = np.array(list(map(lambda x: [len(x)], Df["hit_seq"])))
        return np.absolute(hit_seq - query_seq), hit_seq, query_seq

    def net_charge_calc(self, sequence: pd.Series, pH: int) -> np.array:
        PKA_values = pd.DataFrame(
            data=[
                [2.34, 2.34, 2.32, 2.36, 2.36, 2.28, 1.99, 1.83, 2.83, 2.02, 2.17, 2.21, 2.09, 2.2, 1.96, 1.88, 2.19, 2.18, 2.17, 1.82],
                [9.6, 9.69, 9.62, 9.6, 9.6, 9.21, 10.6, 9.13, 9.39, 8.8, 9.13, 9.15, 9.1, 9.11, 8.18, 9.6, 9.67, 8.95, 9.04, 9.17],
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                 np.nan, np.nan, np.nan, 10.07, 8.33, 3.65, 4.25, 10.53, 12.48, 6]
            ],
            columns= ["G", "A", "V", "L", "I", "M", "P", "F", "W", "N",
                      "Q", "S", "T", "Y", "C", "D", "E", "K", "R", "H"]
        )
        # Henderson_Haselbach
        find_N_terminal = lambda x: PKA_values.loc[1, x[0]]
        find_c_terminal = lambda x: PKA_values.loc[0, x[-1]]
        find_alifate = lambda x: PKA_values.loc[2, x]
        Hend_base = lambda y,x: y.count(x) * (
                    np.power(10, find_alifate(x)) / (np.power(10, pH) + np.power(10, find_alifate(x))))
        Hend_acid = lambda y,x: y.count(x) * (
                    np.power(10, pH) / (np.power(10, pH) + np.power(10, find_alifate(x))))
        Hend_C = lambda x: (10 ** pH) / (10 ** pH + 10 ** find_c_terminal(x))
        Hend_N = lambda x: np.power(10, find_N_terminal(x)) / (np.power(10, pH) + np.power(10, find_N_terminal(x)))
        seq_val = pd.DataFrame(
            [list(map(Hend_N, sequence)),
             list(map(Hend_C, sequence)),
             list(map(Hend_acid, sequence, ['D' for i in range(len(sequence))])),
             list(map(Hend_acid, sequence, ['E' for i in range(len(sequence))])),
             list(map(Hend_acid, sequence, ['C' for i in range(len(sequence))])),
             list(map(Hend_acid, sequence, ['Y' for i in range(len(sequence))])),
             list(map(Hend_base, sequence, ['R' for i in range(len(sequence))])),
             list(map(Hend_base, sequence, ['K' for i in range(len(sequence))])),
             list(map(Hend_base, sequence, ['H' for i in range(len(sequence))]))
            ]
        )
        calculate = np.subtract(np.sum(seq_val.loc[[0,6,7,8],:]), np.sum(seq_val.loc[[1,2,3,4,5],:]))
        return calculate

    def iso_elec_point(self, *args):
        sequence = args[0]
        if len(args) >= 2:
            pH = args[1]
        else:
            pH = np.arange(0, 14, 0.01, dtype=np.float64)

        ph_content = {}
        ph_calc = {}
        for item in pH:
            ph_content[item] = self.net_charge_calc(sequence, item)
        ph_content = pd.DataFrame(ph_content)
        ph_dataframe = np.absolute(ph_content)
        for seq in ph_dataframe.T:
            seq_in_ph = ph_dataframe.loc[seq, :]
            ph_calc[seq] = np.mean(seq_in_ph[seq_in_ph < 0.5].index)
        return pd.Series(ph_calc)

    def sort_by_anno(self, file_name):
        anno = {}
        for probe in self.probes:
            df = pd.read_excel(self.find_files(probe, "{}".format(file_name)), engine='odf', index_col=0)
            anno_ids = tuple(set(df["hit_id"]))
            addon = lambda x: df.query("hit_id == @x")
            sorted_addons = list(map(addon, anno_ids))
            for indx, data_frame in enumerate(sorted_addons):
                data_frame.loc[:, "probe"] = [probe for item in range(len(data_frame))]
                if 0 == self.probes.index(probe):
                    anno[anno_ids[indx]] = data_frame
                else:
                    anno[anno_ids[indx]] = pd.concat([ anno[anno_ids[indx]] , data_frame], ignore_index=True)
        return anno

    def generate_fasta(self, sorted_file_loc, fasta_location):
        df = pd.read_excel(sorted_file_loc)
        names_funk = lambda x: f">{df.query_id[x]} ({df.probe[x]}) \n"
        names = list(map(names_funk, range(len(df))))
        names.append(">{} {} \n".format(df.hit_id[0], df.hit_description[0]))
        sequences = []
        for indx, seq in enumerate(df["query_seq"]):
            i = 0
            seq_str = ""
            while i < len(df["query_seq"][indx]):
                j = i + 80
                seq_str = seq_str + f"{df.query_seq[indx][i:j]}\n"
                i = i + 80
            sequences.append(seq_str)
        i = 0
        hit_seq = ""
        while i < len(df["hit_seq"][0]):
            j = i + 80
            hit_seq = hit_seq+ f"{df.hit_seq[0][i:j]}\n"
            i = i + 80
        sequences.append(hit_seq)
        print(len(sequences), len(names))
        with open(fasta_location, "w") as my_file:
            for indx in range(len(names)):
                my_file.write(names[indx])
                my_file.write(sequences[indx])



    def table_of_occurance(self):
        pass

    def best_sequences_by_HMM(self, df, e_value_tresh): #smaller evalue better
        return df[df['e-value'] < e_value_tresh]

    def cut_to_long_seqs(self, df, max_length): # end of histogram
        new_df = df[df['query_seq_len'] <= max_length]
        return new_df

    def most_active(self):
        #sorting_data.py
        #delete annotations with no occurances
        pass

    def histograms(self, df, df_name):
        plt.hist(df['query_seq_len'], bins=90)
        plt.title("{}".format(df_name))
        plt.xlabel("Sequence lengths")
        plt.ylabel("Number of sequences")
        plt.axvline(x=380, color="red", linestyle='--', label="Cutting point{}".format(380))
        plt.show()

    def RMSD(self, location, traj_nc, traj_prmtop):
        traj = md.load(os.path.join(location, traj_nc), top=os.path.join(location, traj_prmtop))
        prot_sel = traj.topology.select('protein')
        prot_traj = traj.atom_slice(prot_sel)

        rmsds = md.rmsd(prot_traj, prot_traj, 0)
        plt.plot(np.arange(0, len(rmsds)), rmsds, color='m', label='Rep 1')
        plt.title("RMSD of {}".format("Test sequence"))
        plt.ylabel("Å")
        plt.show()
        return rmsds

    def RMSF(self, location, traj_nc, traj_prmtop):
        traj = md.load(os.path.join(location, traj_nc), top=os.path.join(location, traj_prmtop))
        prot_sel = traj.topology.select('protein')
        prot_traj = traj.atom_slice(prot_sel)
        rmsf = md.rmsf(prot_traj, prot_traj,0)
        plt.plot(np.arange(0, len(rmsf)), rmsf, color='g', label='Rep 1')
        plt.title("RMSF of {}".format("Test sequence"))
        plt.ylabel("Å")
        plt.xlabel("Aminoacid")
        plt.show()
        return rmsf
