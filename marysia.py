import pandas
import numpy

#sprawdza czy wyniki surowe piotra sa maja takie same dopasowania czy inne
wyniki = pandas.read_excel("/home/arjissuan/Desktop/AAT/out.ods", engine="odf")
piotr_wyniki = pandas.read_csv("/home/arjissuan/Desktop/AAT/metrics/D7_metrics.csv")

qtitle = piotr_wyniki["qtitle"].str.split(" ")
title = piotr_wyniki["title"].str.split(" ")
print(qtitle[0][0],
      title[0][1])

qtitle = pandas.Series(qtitle[item][0] for item in range(len(qtitle)))
title = pandas.Series(title[ite][1] for ite in range(len(title)))
print(qtitle[0], title[0])

zip_wyniki = pandas.Series((a,b) for a, b in zip(wyniki["query_id"], wyniki["hit_id"]))
piotr_zip_w = pandas.Series((c,d) for c, d in zip(qtitle, title))
print(zip_wyniki[1], piotr_zip_w[99])

wyniki["is_same"] = zip_wyniki.isin(piotr_zip_w)
wyniki["is_same"] = wyniki["is_same"].astype(str)
wyniki["is_same"] = wyniki["is_same"].replace("True", "1")
wyniki["is_same"] = wyniki["is_same"].replace("False", "0")
wyniki["is_same"] = wyniki["is_same"].astype(int)
wyniki.head()
print(type(wyniki["is_same"][1]))
wyniki.to_excel("/home/arjissuan/Desktop/D7/out.ods")



