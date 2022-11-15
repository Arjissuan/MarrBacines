from Bio import SeqIO
import pandas
import numpy

def metrik():
    record = list(SeqIO.parse("/home/arjissuan/Desktop/SD1/mono_out.fsa", "fasta"))
    record = pandas.Series(record[item].id for item in range(len(record)))
    print(record)
    SD1_record = pandas.read_excel("/home/arjissuan/Desktop/SD1/out.ods")
    DBL_record = pandas.read_excel("/home/arjissuan/Desktop/DBL/out.ods")
    D7_record = pandas.read_excel("/home/arjissuan/Desktop/D7/out.ods")
    MSI_record = pandas.read_excel("/home/arjissuan/Desktop/MSI/out.ods")
    MKC_record = pandas.read_excel("/home/arjissuan/Desktop/MKC/out.ods")
    MD89A_record = pandas.read_excel("/home/arjissuan/Desktop/MD89A/out.ods")
    MD87_record = pandas.read_excel("/home/arjissuan/Desktop/MD87/out.ods")
    IK1_record = pandas.read_excel("/home/arjissuan/Desktop/IK1/out.ods")
    M11_record = pandas.read_excel("/home/arjissuan/Desktop/M11/out.ods")
    all_records = pandas.DataFrame(
        columns=["baza_danych", "SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11"])

    all_records["baza_danych"] = pandas.concat([all_records["baza_danych"], record], ignore_index=True)

    for item in range(len(all_records["baza_danych"])):
        value = all_records["baza_danych"][item]
        powt_SD1 = len(SD1_record.query("hit_id == @value"))
        powt_dbl = len(DBL_record.query("hit_id == @value"))
        powt_d7 = len(D7_record.query("hit_id == @value"))
        powt_MSI = len(MSI_record.query("hit_id == @value"))
        powt_MKC = len(MKC_record.query("hit_id == @value"))
        powt_MD89A = len(MD89A_record.query("hit_id == @value"))
        powt_md87 = len(MD87_record.query("hit_id == @value"))
        powt_ik1 = len(IK1_record.query("hit_id == @value"))
        powt_m11 = len(M11_record.query("hit_id == @value"))
        # print(value)
        # print(len(SD1_record.query("hit_id == @value")))
        all_records["SD1"][item] = powt_SD1
        all_records["DBL"][item] = powt_dbl
        all_records["D7"][item] = powt_d7
        all_records["MSI"][item] = powt_MSI
        all_records["MKC"][item] = powt_MKC
        all_records["MD89A"][item] = powt_MD89A
        all_records["MD87"][item] = powt_md87
        all_records["IK1"][item] = powt_ik1
        all_records["M11"][item] = powt_m11
    # print(all_records)
    # all_records.to_csv("/home/arjissuan/Desktop/tavle/table_with_numbers.csv")

#tworzenie metryki ilosci wystapien annotacji z prokka(hitow) w wynikach
def aktywnosc(lista_prob):
    anotations = list(SeqIO.parse("/home/arjissuan/Desktop/SD1/mono_out.fsa", "fasta"))
    anotations = pandas.Series(anotations[item].id for item in range(len(anotations)))
    metric = pandas.DataFrame()
    metric["baza_danych"] = anotations
    for proba in lista_prob:
        proba_records = pandas.read_excel(r"/home/arjissuan/Desktop/{}/out_with_low_eV.ods".format(proba), engine='odf')
        print(proba_records)
        proba_sumy = []
        for item in range(len(metric["baza_danych"])):
            annotation_id = metric["baza_danych"][item]
            suma_powt = len(proba_records.query("hit_id == @annotation_id"))
            proba_sumy.append(suma_powt)
        metric[proba] = proba_sumy
        print(metric)

    return metric


metr = aktywnosc(["SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11", "MD19A", "AAT"])
metr.to_excel(r"/home/arjissuan/Desktop/tavle/tablica_wysta.ods")



