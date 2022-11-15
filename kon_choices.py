from Bio import SearchIO, SeqIO
import pandas
#2
#tworzenie parsowalnej tabeli wynikow out.ods dla kazdej proby
file1 = "/home/arjissuan/Desktop/MD19A/results_MD19A.txt"
dane = []
for item in SearchIO.parse(file1, "phmmer3-domtab"):
    dane.append(item)
#works only with phmmer --domtab, SearchIO
#dane[0] query
#dane[0][0] hit
#dane[0][0][0] hsp?
#dane[0][0][0][0] hspfragment
#print(dane[0][0].query_id)
#print(dane[0][0].id)
#print(dane[0][0].bitscore)

#making working database of sequences
record_que = list(SeqIO.parse("/home/arjissuan/Desktop/MD19A/PROKKA_all.faa","fasta"))
record_que_dict = {"name":[], "seq":[]}
for item in range(len(record_que)):
    record_que_dict["name"].append(record_que[item].name)
    record_que_dict["seq"].append(record_que[item].seq)

record_hit = list(SeqIO.parse("/home/arjissuan/Desktop/MD19A/mono_out.fsa", "fasta"))
record_hit_dict = {"name":[], "seq":[]}
for item in range(len(record_hit)):
    record_hit_dict["name"].append(record_hit[item].name)
    record_hit_dict["seq"].append(record_hit[item].seq)

#data from results to dictionary
data_readable = {"query_id":[], "query_seq":[], "bitscore":[],'bias':[], "e_value":[], "hit_id":[], "hit_seq":[]}
for ite in range(len(dane)):
    for i in range(len(dane[ite])):
        print(dane[ite][i])
        data_readable["query_id"].append(dane[ite][i].query_id)
        data_readable["bitscore"].append(dane[ite][i].bitscore)
        data_readable["e_value"].append(dane[ite][i].evalue)
        data_readable["hit_id"].append(dane[ite][i].id)
        data_readable["bias"].append(dane[ite][i].bias)

#looking for sequences of querries and hits
for item in range(len(data_readable["query_id"])):
    if data_readable["query_id"][item] in record_que_dict["name"]:
        indx = record_que_dict["name"].index(data_readable["query_id"][item])
        print(indx)
        data_readable["query_seq"].append(record_que_dict["seq"][indx])

for item in range(len(data_readable["hit_id"])):
    if data_readable["hit_id"][item] in record_hit_dict["name"]:
        indxx = record_hit_dict["name"].index(data_readable["hit_id"][item])
        print(indxx)
        data_readable["hit_seq"].append(record_hit_dict["seq"][indxx])





print(len(data_readable["query_id"]),
      len(data_readable["query_seq"]),
      len(data_readable["hit_id"]),
      len(data_readable["hit_seq"]),)

dt = pandas.DataFrame(data_readable)
dt.to_excel("/home/arjissuan/Desktop/MD19A/out.ods")















