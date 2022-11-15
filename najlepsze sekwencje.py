import pandas as pd

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

Threshold().filtrowanie_po_anotacji(probki)



