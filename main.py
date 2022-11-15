from Bio import SeqIO
#1
#odzielanie bialek o wielu sekwencjach od tych z jedna
record = list(SeqIO.parse("/home/arjissuan/Desktop/SD1/out.fsa", "fasta"))
print(record[0].seq)
if ',' not in record[0].seq:
    print(True)


def cleaning(record):
    mono = []
    poli = []
    for item in range(len(record)):
        if "," in record[item].seq:
            poli.append(record[item])
        elif "-" in record[item].seq:
            poli.append(record[item])
        else:
            mono.append(record[item])

    with open("/home/arjissuan/Desktop/SD1/mono_out.fsa", "w") as file:
        SeqIO.write(mono, file, "fasta")

    with open("/home/arjissuan/Desktop/SD1/poli_out.fsa", "w") as file:
        SeqIO.write(poli, file, "fasta")










