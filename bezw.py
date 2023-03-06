import pandas
import numpy

probki = ["SD1", "DBL", "D7", "MSI", "MKC", "MD89A", "MD87", "IK1", "M11", "MD19A", "AAT"]
annotats = pandas.read_excel("/home/arjissuan/Desktop/tavle_significant/interesujace_anotacje.ods", engine='odf') #['baza_danych']
#obliczania wart bezwlglednej
def bezwzgledna(frame):
    df = pandas.read_excel(frame)
    qlen = df["query_seq"]  # type pandas.series from dataframe
    hlen = df["hit_seq"]

    wart = numpy.array(qlen.str.len()) - numpy.array(hlen.str.len())
    seria = pandas.Series(numpy.absolute(wart))
    df["absolute"] = seria
    df.head()
    return df.to_excel(frame)
#bezwzgledna("/home/arjissuan/Desktop/MD19A/out.ods")


def del_of_longer(frame):
    df = pandas.read_excel(frame, engine='odf')
    df = df.drop(columns=["Unnamed: 0"])
    qlen = df['query_seq']
    hlen = len(df['hit_seq'][0])
    tru_len = int(hlen*3)
    print(tru_len, hlen, hlen*1.5)

    dframe = (df[qlen.map(lambda x: len(str(x))) < 100])
    if tru_len > 100:
        newer_ver = (dframe[qlen.map(lambda x: len(str(x))) < hlen])
    elif tru_len <= 100:
        newer_ver = (dframe[qlen.map(lambda x: len(str(x))) < tru_len])
    else:
        newer_ver = pandas.DataFrame(columns=list(df.columns))

    return newer_ver




if __name__ == '__main__':
    for item in annotats['baza_danych']:
        df1 = del_of_longer("/home/arjissuan/Desktop/tavle_significant/best_seq_for_anno_{}.ods".format(item))
        df1.to_excel("/home/arjissuan/Desktop/tavle_significant/after_deltion_of_too_long_seq/{}_after_erasing_long_seq.ods".format(item))
