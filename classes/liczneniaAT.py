from Bio import SeqIO
import numpy
import pandas
a = "out.fsa"
b = "/home/arjissuan/Downloads/myco.fasta"


def at_counter(path):
    sekw = SeqIO.read(r"{}".format(path), 'fasta')
    arr_sekw = pandas.Series(data=list(sekw.seq), name="Sequences")
    arr_sekw = arr_sekw.replace(["A", "T", "C", "G"], [11, 2, 3, 9])
    trve_sekw = arr_sekw.value_counts()
    indexy = numpy.array(pandas.Series(trve_sekw.index).replace([11, 2, 3, 9], ["A", "T", "C", "G"]))
    seria = pandas.Series(data=trve_sekw.values, index=indexy)
    try:
        AA = numpy.divide((seria["A"] * 100), seria.sum())
    except KeyError:
        AA = 0
    try:
        TT = numpy.divide((seria["T"] * 100), seria.sum())
    except KeyError:
        TT = 0
    return ("Zawartość każdego nukleotydu:\n{}\nZawartosc A: {}\nZawartosc T: {}\nSuma zawartosci procentowej{}\nSuma AT {}".format(seria, AA, TT, AA+TT, seria["T"]+seria["A"]))



print(at_counter(b))