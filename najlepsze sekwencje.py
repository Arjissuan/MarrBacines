import pandas as pd
from Bio import SearchIO, SeqIO

probki = ["SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11", "MD19A", "AAT"]
class Threshold:
    plik = pd.read_excel(r"/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods", engine='odf')
    def filtrowanie_po_probie(self, lista_probek):
        annotations = self.plik["baza_danych"]
        for prob in lista_probek:
            data = pd.read_excel(r"/home/arjissuan/Desktop/{}/out_with_low_eV.ods".format(prob), engine='odf')
            most_active = pd.DataFrame()
            for id in annotations:
                sig_data = data.query("hit_id == @id")
                most_active = pd.concat([most_active, sig_data])
            most_active = most_active.drop(columns="Unnamed: 0")
            print(most_active)
            most_active.to_excel(r"/home/arjissuan/Desktop/{}/Most_active_from_{}.ods".format(prob, prob))

    def filtrowanie_po_anotacji(self, lista_probek):
        annotations = self.plik["baza_danych"]
        for anotat in annotations:
            zbiornik = pd.DataFrame()
            for probe in lista_probek:
                data = pd.read_excel("/home/arjissuan/Desktop/{}/out_with_low_eV.ods".format(probe) , engine='odf')
                fetch = data.query("hit_id == @anotat")
                klucz = pd.DataFrame()
                if len(fetch)>0:
                    klucz = pd.concat([klucz, fetch])
                    id_probe = []
                    for item in range(len(fetch)):
                        id_probe.append(probe)
                    klucz["id_probe"] = id_probe
                    zbiornik = pd.concat([zbiornik, klucz])
                    zbiornik = zbiornik.drop(columns="Unnamed: 0")
            print(zbiornik)
            zbiornik.to_excel(r"/home/arjissuan/Desktop/tavle_significant/best_seq_for_anno_{}.ods".format(anotat))

    def to_FASTA(self,probes):
        ln = 60
        data = pd.read_excel("/home/arjissuan/Desktop/tavle_significant/best_seq_for_anno_{}.ods".format(probes), engine='odf')
        q_id = data["query_id"]
        q_seq = data["query_seq"]
        query = zip(q_id, q_seq)
        fasta_whole = ""
        for item in query:
            sekw = item[1]
            if len(sekw) <= 60:
                fasta_rec = ">{}\n{}\n".format(item[0], sekw)
                fasta_whole+= fasta_rec
            elif len(sekw) > 60:
                sekw_zb = ""
                while len(sekw) !=0:
                    sekw_zb = sekw_zb + "{}\n".format(sekw[:60])
                    sekw = sekw[60:]
                fasta_re = ">{}\n{}".format(item[0], sekw_zb)
                fasta_whole+=fasta_re
        return fasta_whole




#Threshold().filtrowanie_po_anotacji(probki)
#Threshold().filtrowanie_po_probie(probki)
for item in pd.read_excel("/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods")["baza_danych"]:
    with open("/home/arjissuan/PycharmProjects/pythonProject/annotations_fastas/bestseqfor_{}.fasta".format(item), "w") as file:
        file.write(Threshold().to_FASTA(item))



