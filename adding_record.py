import pandas
from Bio import SeqIO, SearchIO
#do dodawania rzeczy do istniejacych juz tabel aktywnosci

class Adding:
    probki = ["SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11", "MD19A", "AAT"]
    def add(self, rekord, sec_path):
        path = "/home/arjissuan/Desktop/tavle/table_with_numbers.ods"
        base = pandas.read_excel(path)
        read_probe = pandas.read_excel(sec_path , engine='odf')
        ilosc = []
        for item in range(len(base)):
            value = base["baza_danych"][item]
            powt = len(read_probe.query("hit_id == @value"))
            ilosc.append(powt)
        base[rekord] = ilosc
        base.to_excel(path)
        return "Done"
#dodawanie tabelek do wynikow czyli out.ods
    def add_bias_to_out(self, txt, ods):
        dane = []
        for item in SearchIO.parse(txt, "phmmer3-domtab"):
            dane.append(item)
        kolu = []
        for item in range(len(dane)):
            for dex in range(len(dane[item])):
                kolu.append(dane[item][dex].bias)
        file = pandas.read_excel(ods, engine="odf")
        file["bias"] = kolu
        return file

    def add_descriptions_to_existing(self):
        plik = list(SeqIO.parse("/home/arjissuan/Desktop/SD1/mono_out.fsa", "fasta"))
        plik_des = {"Deskrypcja":[], "hit_id":[]}
        for item in plik:
            #print(item.description)
            plik_des["Deskrypcja"].append(item.description)
            plik_des["hit_id"].append(item.name)
        plik_de_pd = pandas.DataFrame(plik_des)
        print(plik_de_pd)
        tablica_inter_annota = pandas.read_excel("/home/arjissuan/Desktop/tavle_significant/tablica_wysta.ods", engine="odf")
        bufor = pandas.DataFrame()
        for item in range(len(tablica_inter_annota)):
            value = tablica_inter_annota["baza_danych"][item]
            a = plik_de_pd.query("hit_id == @value")
            bufor = pandas.concat([bufor, a])
        new = []
        for item in bufor["Deskrypcja"]:
            new.append(item)
        tablica_inter_annota["Descritpion"] = new
        print(tablica_inter_annota)
        return tablica_inter_annota

    def add_decription_to_aligments(self, probe):
        file = pandas.read_excel(r"/home/arjissuan/Desktop/{}/out_with_low_eV.ods".format(probe,probe), engine="odf")
        file = file.drop(columns="Unnamed: 0")
        des = list(SeqIO.parse(r"/home/arjissuan/Desktop/{}/PROKKA_all.faa".format(probe), 'fasta'))
        des_list = list((item.description, item.name) for item in des)
        bufor = pandas.DataFrame()
        decri = []
        for dest_item in des_list:
            value = dest_item[1]
            if len(file.query("query_id == @value")) > 0:
                bufor = pandas.concat([bufor, file.query("query_id == @value")])
                for i in range(len(file.query("query_id == @value"))):
                    decri.append(dest_item[0])

        bufor["query_description"] = decri
        print(bufor)
        return bufor



#Adding().add("D7", "/home/arjissuan/Desktop/D7/out.ods")
#print(pandas.read_excel("/home/arjissuan/Desktop/AAT/out.ods", engine='odf'))

#new_file = (Adding().add_bias_to_out("/home/arjissuan/Desktop/AAT/results_AAT.txt", "/home/arjissuan/Desktop/AAT/out.ods"))
#new_file.to_excel("/home/arjissuan/Desktop/AAT/out.ods")

#Adding().add_descriptions_to_existing().to_excel("/home/arjissuan/Desktop/tavle_significant/tablica_wysta_update.ods")
for probe in Adding.probki:
    Adding().add_decription_to_aligments(probe).to_excel(r"/home/arjissuan/Desktop/{}/out_with_low_eV.ods".format(probe))
