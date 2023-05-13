import numpy
import pandas
#usuwa wyniki ze zbyt wysokim evalue z wynikow /{}/out.ods
#table = pandas.read_excel("/home/arjissuan/Desktop/tavle/table_with_numbers.ods", engine='odf')
def numerous(list_of_probes):
    for probe in list_of_probes:
        tb = pandas.read_excel(r'/home/arjissuan/Desktop/{}/out.ods'.format(probe), engine='odf')
        arr = tb["e_value"]
        cond = []
        for i in range(len(arr)):
            cond.append(0.0003)

        after = numpy.less_equal(arr, cond)
        tb["condition"] = after.astype(str)
        tb["condition"] = tb["condition"].replace("True", "1")
        tb["condition"] = tb["condition"].replace("False", "0")
        tb["condition"] = tb["condition"].astype(int)
        tb = tb.query("condition == 1")
        tb.to_excel(r'/home/arjissuan/Desktop/{}/out_with_low_eV.ods'.format(probe))

#numerous(["SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11", "MD19A", "AAT"])

#funkcja ktora posortuje metryczke oraz usunie annotacje zerowe(bez aktywnosci)
def significant(matrix):
    tavle = pandas.read_excel(r"/home/arjissuan/Desktop/tavle/{}".format(matrix), engine='odf')
    table = tavle.to_numpy()
    wart = table[:, 2:]
    suma = numpy.array(0)
    for itme in range(len(wart[0, :])):
        print(wart[:, itme])
        suma = suma + wart[:, itme]
    arr = []
    for item in range(len(suma)):
        arr.append(0)

    suma = numpy.greater(suma, arr)
    tavle["suma"] = suma
    sig_tavle = tavle.query("suma == True")
    sig_tavle.to_excel(r"/home/arjissuan/Desktop/tavle_significant/{}".format(matrix))

#significant('table_with_numbers.ods')
#
def most_active(matrix):
    tavle = pandas.read_excel(r"/home/arjissuan/Desktop/tavle_significant/{}".format(matrix), engine='odf')
    table = tavle.to_numpy()
    wart = table[:, 2:]
    suma = numpy.array(0)
    for itme in range(len(wart[0, :])):
        print(wart[:, itme])
        suma = suma + wart[:, itme]

    tavle["suma"] = suma
    info = tavle.sort_values(by="suma")
    info.to_excel(r"/home/arjissuan/Desktop/tavle_significant/{}".format(matrix))

#most_active("tablica_wysta.ods")

for file in ['table_with_numbers.ods', 'tablica_wysta.ods']:
    #significant(file)
    most_active(file)