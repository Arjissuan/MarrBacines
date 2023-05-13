import pandas
import pandas as pd
#skrypt tworzacy tabele dla surowych wynikow, ktora zawiera wylacznie rozne wyniki(hity, out.fsa) dla tych samych zapytan(query(anotacji z prokka))
wyniki = pd.read_excel("/home/arjissuan/Desktop/M11/out.ods", engine="odf")
nie_pas_dop = pd.DataFrame({})
q_id = []
q_seq = []
h_id = []
h_seq = []
for item in range(len(wyniki)):
    if wyniki["is_same"][item] == 0:
        q_id.append(wyniki["query_id"][item])
        q_seq.append(wyniki["query_seq"][item])
        h_id.append(wyniki["hit_id"][item])
        h_seq.append(wyniki["hit_seq"][item])

nie_pas_dop["query_id"] = q_id
nie_pas_dop["query_seq"] = q_seq
nie_pas_dop["hit_id"] = h_id
nie_pas_dop["hit_seq"] = h_seq
print(nie_pas_dop)

#tworze DT dla tych samych zapoytan ale roznych query w odniesieniu do wynikow z blasta
piotr_wyniki = pd.read_csv("/home/arjissuan/Desktop/M11/metrics/M11_metrics.csv")
print(piotr_wyniki.columns)
qtitle = piotr_wyniki["qtitle"].str.split(" ")
title = piotr_wyniki["title"].str.split(" ")
title_seq = piotr_wyniki["hitseq"]

qtitle = pd.Series(qtitle[item][0] for item in range(len(qtitle)))
title = pd.Series(title[item][1] for item in range(len(title)))

z_p_w_q = []
z_p_w_h = []
z_p_w_hs = []

Q_bools = qtitle.isin(nie_pas_dop["query_id"])
H_bools = title.isin(nie_pas_dop["hit_id"])
HQ_bools = pandas.DataFrame({"Q_bools":Q_bools, "H_bools":H_bools})
print(HQ_bools)
for item in range(len(HQ_bools)):
    if HQ_bools["Q_bools"][item] == True and HQ_bools["H_bools"][item] == False:
        z_p_w_q.append(qtitle[item])
        z_p_w_h.append(title[item])
        z_p_w_hs.append(title_seq[item])

z_p_w = pandas.DataFrame({"query_id":z_p_w_q, "hit_id_blast":z_p_w_h, "hit_seq_blast":z_p_w_hs})#tworze zbior interesujacych wynikow piotra
print(z_p_w)
set_of_query = set(q_id)#unikalne id
dane = pandas.DataFrame({})
for item in set_of_query:
    value = item
    moje = nie_pas_dop.query("query_id == @value")
    dane = pandas.concat([dane, moje], ignore_index=True)
    piotra = z_p_w.query("query_id == @value")
    dane = pd.concat([dane, piotra], ignore_index=True)


print(dane)
print(dane.columns)
dane.to_csv("/home/arjissuan/Desktop/M11/metrykaM11.csv")


