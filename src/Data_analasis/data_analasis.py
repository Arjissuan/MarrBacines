from Bio import SeqIO, SearchIO
import pandas as pd
import numpy as np
import os


class MarrBacines:
    def __init__(self, location="/home/arjissuan/Desktop/", probes=("DBL", "IK1", "M11", "MD19A", "MD87", "MD89A", "MKC", "MSI", "SD1", "D7", "AAT")):
        self.probes = probes
        self.location = location
        self.DB_name = "PROKKA_all.faa"

    def DB_prep_for_phmmer(self):
        #probably main
        pass

    def find_files(self, probe, results_name):
        files = os.listdir(self.location + probe)
        name_of_res = lambda x: x if x in results_name else None
        name = str(list(map(name_of_res, files)).remove(None))
        return name

    def domtab_to_excel(self, *results_name):
        #kon_choices
        for probe in self.probes:
            if results_name == ():
                file = self.location + probe + f"/results_{probe}.txt"
            else:
                file = self.find_files(probe, results_name[0])
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

            funk = lambda x: (x.name, x.seq)
            AMP_DB = dict(map(funk, SeqIO.parse(self.location+probe+"/mono_out.fsa", "fasta")))
            AMP_probe = dict(map(funk, SeqIO.parse(self.location+probe+"/PROKKA_all.faa", "fasta")))

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

            df.to_excel(self.location+probe+f"/out_{probe}.ods")

