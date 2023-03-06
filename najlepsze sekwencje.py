import pandas as pd



class Threshold:
    def __init__(self, lista_probek):
        self.lista_probek = lista_probek

    plik = pd.read_excel(r"/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods", engine='odf')
    def filtrowanie_po_probie(self):
        annotations = self.plik["baza_danych"]
        for prob in self.lista_probek:
            data = pd.read_excel(r"/home/arjissuan/Desktop/{}/out_with_low_eV.ods".format(prob), engine='odf')
            most_active = pd.DataFrame()
            for id in annotations:
                sig_data = data.query("hit_id == @id")
                most_active = pd.concat([most_active, sig_data])
            most_active = most_active.drop(columns="Unnamed: 0")
            print(most_active)
            most_active.to_excel(r"/home/arjissuan/Desktop/{}/Most_active_from_{}.ods".format(prob, prob))

    def filtrowanie_po_anotacji(self):
        annotations = self.plik["baza_danych"]
        for anotat in annotations:
            zbiornik = pd.DataFrame()
            for probe in self.lista_probek:
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

    #makes fasta files from sequence aligned to each annoted(from DB)
    def to_FASTA(self,probes):
        ln = 60
        data = pd.read_excel("/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/{}_after_erasing_long_seq.ods".format(probes), engine='odf')
        q_id = data["query_id"]
        q_seq = data["query_seq"]
        q_probe = data["id_probe"]
        query = zip(q_id, q_seq, q_probe)

        fasta_whole = ""
        for item in query:
            sekw = item[1]
            if len(sekw) <= 60:
                fasta_rec = ">{} ; {}\n{}\n".format(item[0], item[2], sekw)
                fasta_whole+= fasta_rec
            elif len(sekw) > 60:
                sekw_zb = ""
                while len(sekw) != 0:
                    sekw_zb = sekw_zb + "{}\n".format(sekw[:60])
                    sekw = sekw[60:]
                fasta_re = ">{} ; {}\n{}\n".format(item[0], item[2], sekw_zb)
                fasta_whole += fasta_re

        h_id = data['hit_id'][0]
        h_seq = data["hit_seq"][0]

        if len(h_seq) <= 60:
            fasta_part = ">{}\n{}\n".format(h_id, h_seq)
            fasta_whole += fasta_part
        elif len(h_seq) > 60:
            sekw = h_seq
            sekw_zb = ""
            while len(sekw) != 0:
                sekw_zb = sekw_zb + "{}\n".format(sekw[:60])
                sekw = sekw[60:]
            fasta_re = ">{}\n{}".format(h_id, sekw_zb)
            fasta_whole += fasta_re

        return fasta_whole

    #makes fasta of probes sequence in relation to annoted sequences. Which means they will be one probe to one annotation sequence
    def to_Fasta_for_probes(self, probe):
        path = r"/home/arjissuan/Desktop/{}/Most_active_from_{}.ods".format(probe, probe)
        data = pd.read_excel(path, engine='odf')
        for anno in pd.read_excel("/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods")["baza_danych"]:
            fasta_whole = ""
            value = anno
            carnister = data.query("hit_id == @value")
            if len(carnister) > 0:
                zipped = zip(carnister["query_id"], carnister["query_seq"])
                for item in zipped:
                    if len(item[1]) <= 60:
                        fasta = ">{}\n{}\n".format(item[0], item[1])
                    elif len(item[1]) > 60:
                        seq = item[1]
                        sekw = ""
                        while len(seq) > 0:
                            sekw +="{}\n".format(seq[:60])
                            seq = seq[60:]
                        fasta = ">{}\n{}".format(item[0], sekw)
                    fasta_whole+=fasta
                with open("/home/arjissuan/Desktop/{}/probe_seq_for_anno/seq_for_{}.fasta".format(probe, value), "w") as file:
                    file.write(fasta_whole)
                print(fasta_whole)
            else:
                fasta_whole = f"Dla {value}w probie {probe} nie ma sekwencji"
                with open("/home/arjissuan/Desktop/{}/probe_seq_for_anno/seq_for_{}.txt".format(probe, value), "w") as file:
                    file.write(fasta_whole)
                print(fasta_whole)


probki = ["SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11", "MD19A", "AAT"]

#Pliki = Threshold(probki)
#Pliki.filtrowanie_po_anotacji()
#Pliki.filtrowanie_po_probie()
#for item in pd.read_excel("/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods")["baza_danych"]:
#    with open("/home/arjissuan/PycharmProjects/Pandas+numpy/annotations_fastas/bestseqfor_{}.fasta".format(item), "w") as file:
#        file.write(Threshold().to_FASTA(item))

#for item in probki:
#    Threshold().to_Fasta_for_probes(item)




annotats = pd.read_excel("/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods", engine='odf')

if __name__ == '__main__':
    for item in annotats['baza_danych']:
        with open("/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/{}_a_e_l_s.fasta".format(item),
                  "w") as my_file:
            my_file.write(Threshold(probki).to_FASTA(item))