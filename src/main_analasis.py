from Data_analasis.data_analasis import MarrBacines
import pandas as pd
import os
import numpy as np
import time
class Done_tasks:
    def __init__(self, adres):
        self.analizator = MarrBacines(location=adres)

    def from_dom_to_exc(self):
        for names in [("schmidt", "Schmidt17"), ("branson", "branson23")]:
            self.analizator.domtab_to_excel(names[0], names[1])

    def adding_col(self):
        for names in ["schmidt.ods", "branson.ods"]:
            absulute = list(map(lambda x: self.analizator.absolute(x, names), self.analizator.probes))
            new_values = list(map(lambda x: pd.read_excel(self.analizator.find_files(x, names), engine='odf').drop(columns=["Unnamed: 0"]), self.analizator.probes))
            for item in range(len(new_values)):
                new_values[item]["Absolute"] = absulute[item][0]
                new_values[item]["hit_seq_len"] = absulute[item][1]
                new_values[item]["query_seq_len"] = absulute[item][2]
                new_values[item].to_excel(self.analizator.find_files(self.analizator.probes[item], names), engine="odf")

    def chem_prope(self, names=("branson.ods", "schmidt.ods")):
        for name in names:
            list_df = list(map(lambda x: self.analizator.AMP_chem_prop(x, name), self.analizator.probes))
            for indx, item in enumerate(list_df):
                probe = self.analizator.probes[indx]
                item.to_excel(f"/home/arjissuan/Desktop/chemical_properties/{probe}/chem_prop_{name}", engine="odf")

    def sorting(self, *from_sort_by_anno):
        #data from sort by anno
        if 0 == len(from_sort_by_anno):
            file = str(input("branson.ods or schmidt.ods?"))
            from_sort_by_anno = self.analizator.sort_by_anno(file)
        whole_data = from_sort_by_anno
        save_in_files = lambda x: whole_data[x].to_excel(os.path.join(self.analizator.location,
                                                                      "by_annotation",
                                                                      x,
                                                                      f"all_of_{x}_test.ods"))
        mkdir = lambda x: os.mkdir(os.path.join(self.analizator.location, "by_annotation", x))
        for key in whole_data.keys():
            try:
                mkdir(key)
            except FileExistsError:
                print("Directory exists")
            save_in_files(key)

    def selection_of_good_sequences(self, annotation, seq_len, evalue):
        location = f"/home/arjissuan/Desktop/cutinasy/by_annotation/{annotation}/all_of_{annotation}.ods"
        df = pd.read_excel(location, engine='odf', index_col=0)
        df = self.analizator.cut_to_long_seqs(df, seq_len)
        df = self.analizator.best_sequences_by_HMM(df, evalue)
        new_location = f"/home/arjissuan/Desktop/cutinasy/by_annotation/{annotation}/sorted_{annotation}.ods"
        df.to_excel(new_location, engine='odf')
        return df


if __name__ == "__main__":
    task = Done_tasks(adres="/home/arjissuan/Desktop/cutinasy")
    analisys = MarrBacines(location="/home/arjissuan/Desktop/cutinasy")
    annotations = [name for name in os.listdir("/home/arjissuan/Desktop/cutinasy/by_annotation")
                   if os.path.isdir(os.path.join("/home/arjissuan/Desktop/cutinasy/by_annotation", name))]


    start = time.time()
    values = []
    for seq in pd.read_excel("/home/arjissuan/Desktop/cutinasy/chemical_properties/SD1/chem_prop_schmidt.ods", index_col=0)["sequence"]:
        values.append(analisys.iso_elec_point(seq))
    print(values, time.time()- start)


