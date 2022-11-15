import pandas
import numpy

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

bezwzgledna("/home/arjissuan/Desktop/MD19A/out.ods")
with open("/home/arjissuan/Desktop/MD19A/out.ods") as file:
    print(file)
a = pandas.Series([1,2,3,4,5,7,5,7,5,5,5,5,1,2])
b = [1,5]
print(a.isin(b))




