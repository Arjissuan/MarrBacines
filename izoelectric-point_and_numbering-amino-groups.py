import pandas as pd
import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os


class Protein_atributes:

    def izoel_point(self):
        #musi sprawdzac pozycje aminikwasow zanim przydzieli wartosci
        anno = (item for item in SeqIO.parse("/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/12690_a_e_l_s_bez.fasta", 'fasta'))
        matrix_of_pka = pd.DataFrame(index=["R", "K", "H", "Y", "C", "E", "D"],
                                     columns=["pKa1", "pKa2", "pKa3", "pI"],
                                     data=[
                                        [2.1, 9, 12.5, 10.8],
                                        [2.1, 9.8, 3.9, 3],
                                        [1.7, 10.4,8.3, 5.0],
                                        [2.2, 9.7,4.3, 3.2],
                                        [1.8, 9.2, 6, 7.8],
                                        [2.2, 9, 10.5, 9.8],
                                        [2.2, 9.1, 10.1, 5.7]]
                                     )



        return matrix_of_pka

    def izooo(self):
        in_put = "Most_active_from_MD19A"
        data = pd.read_excel(
            r"/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/7715_after_erasing_long_seq.ods",
            engine='odf')
        wn = data["query_seq"].to_numpy().astype("str")
        arg = np.array(list(wn.item(indx).count("R") for indx in range(len(wn))))  # r
        lys = np.array(list(wn.item(indx).count("K") for indx in range(len(wn))))  # k
        his = np.array(list(wn.item(indx).count("H") for indx in range(len(wn))))  # h
        tyr = np.array(list(wn.item(indx).count("Y") for indx in range(len(wn))))  # y
        cys = np.array(list(wn.item(indx).count("C") for indx in range(len(wn))))  # c
        glu = np.array(list(wn.item(indx).count("E") for indx in range(len(wn))))  # e
        asp = np.array(list(wn.item(indx).count("D") for indx in range(len(wn))))  # d

        # print(arg, lys, his, tyr, cys, glu, asp)
        # print(arg)
        # print(np.divide(arg, (1 + np.power(10, (10.46 - 9)))))
        pH = 0
        wektor_izooelectrics = np.array([list(0 for item in range(len(wn))), ])
        # wektor = np.array([list(0 for item in range(len(wn))),])
        while pH <= 14:
            CCter = np.divide(1, (1 + np.power(10, (3.65 - pH))))  # cterminal
            CR = np.divide(arg, (1 + np.power(10, (pH - 12.48))))
            CK = np.divide(lys, (1 + np.power(10, (pH - 10.54))))
            CH = np.divide(his, (1 + np.power(10, (pH - 6.04))))
            CY = np.divide(-tyr, (1 + np.power(10, (10.46 - pH))))
            CC = np.divide(-cys, (1 + np.power(10, (8.18 - pH))))
            CE = np.divide(-glu, (1 + np.power(10, (4.07 - pH))))
            CD = np.divide(-asp, (1 + np.power(10, (3.9 - pH))))
            CNH2 = np.divide(1, (1 + np.power(10, (pH - 8.2))))  # ntrminal
            chrge = np.sum(np.array((CD, CE, CC, CY, CH, CK, CR, CCter, CNH2), dtype=object))
            doitisin = np.array(np.isin(np.around(chrge, 1), [0, -0])).astype(float)
            doitisin[doitisin == 1] = pH
            wektor_izooelectrics = np.concatenate((wektor_izooelectrics, [doitisin]), axis=0)
            # wektor = wektor + doitisin
            pH += 0.001
            print(pH)

        seria = []
        for seqid in range(len(wn)):
            arr = np.array(0).astype(float)
            for indx in range(len(wektor_izooelectrics)):
                arr = np.append(arr, wektor_izooelectrics[indx][seqid])
            izo_value = np.divide(np.sum(arr), len(arr[arr > 0]))
            seria.append(izo_value)
            if len(arr[arr > 0]) == 0:
                print("Zero!")
        if len(seria) == len(wn):
            print(True, len(seria), len(wn), seria)
            return seria
        else:
            print(False)

    def countig_statistcs_of_sequences(self, directory):
        # IsoelectricPoint().pi() zwraca iso_point obiekt,a pi() wyciaga z niego pkt_izoo
        anno = []
        for dir in os.listdir('/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq'):
            d = os.path.join('/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq', dir)
            if os.path.isdir(d):
                anno.append(str(dir))

        sec_lam = lambda x: {"sequence":(str(x.seq)),
                             "id":(x.id),
                             'length':(len(x.seq)),
                             "aromaticity":(ProteinAnalysis(str(x.seq)).aromaticity()),
                             "hydrofobicity":(ProteinAnalysis(str(x.seq)).gravy()),
                             "instability_index":(ProteinAnalysis(str(x.seq)).instability_index()), #stable <40 unstable >40
                             "masa":(ProteinAnalysis(str(x.seq)).molecular_weight()),
                             "sec_structure":(ProteinAnalysis(str(x.seq)).secondary_structure_fraction()),
                             "flex":(ProteinAnalysis(str(x.seq)).flexibility()),
                             "izo_punkt":(IsoelectricPoint(x.seq).pi()),
                             }

        files = list(map(lambda x: directory.format(x), anno)) #directory is string with {} to fill the gap with annotaton
        index = 0
        for file in files:
            sec_wektor = pd.DataFrame(list(map(sec_lam, SeqIO.parse(file, 'fasta'))))
            sec_wektor.to_excel(os.path.join("/home/arjissuan/Desktop/tavle_significant/Matrixes_with_statistics", "{}.ods".format(str(anno[index]))))
            index += 1

        return anno


if __name__ == "__main__":
    a = Protein_atributes()
    #print(a.iso_point("/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/{}_a_e_l_s_bez.fasta"))
    # izooo()
    print(a.izoel_point())
