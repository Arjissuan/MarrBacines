import pandas
import pandas as pd
import numpy as np


def izooo():
    in_put = "Most_active_from_MD19A"
    data = pd.read_excel(r"/home/arjissuan/Desktop/MD19A/Most_active_from_MD19A.ods", engine='odf')
    wn = data["query_seq"].to_numpy().astype("str")
    arg = np.array(list(wn.item(indx).count("R") for indx in range(len(wn))))  # r
    lys = np.array(list(wn.item(indx).count("K") for indx in range(len(wn))))  # k
    his = np.array(list(wn.item(indx).count("H") for indx in range(len(wn))))  # h
    tyr = np.array(list(wn.item(indx).count("Y") for indx in range(len(wn))))  # y
    cys = np.array(list(wn.item(indx).count("C") for indx in range(len(wn))))  # c
    glu = np.array(list(wn.item(indx).count("E") for indx in range(len(wn))))  # e
    asp = np.array(list(wn.item(indx).count("D") for indx in range(len(wn))))  # d

    #print(arg, lys, his, tyr, cys, glu, asp)
    #print(arg)
    #print(np.divide(arg, (1 + np.power(10, (10.46 - 9)))))
    pH = 0
    wektor_izooelectrics = np.array([list(0 for item in range(len(wn))),])
    #wektor = np.array([list(0 for item in range(len(wn))),])
    while pH <= 14:
        CCter = np.divide(1, (1 + np.power(10, (3.65 - pH))))  # cterminal
        CR = np.divide(arg, (1 + np.power(10, (pH - 12.48))))
        CK = np.divide(lys, (1 + np.power(10, (pH - 10.54))))
        CH = np.divide(his, (1 + np.power(10, (pH - 6.04))))
        CY = np.divide(-tyr, (1 + np.power(10, (10.46 - pH))))
        CC = np.divide(-cys, (1 + np.power(10, (8.18 - pH))))
        CE = np.divide(-glu, (1 + np.power(10, (4.07 - pH))))
        CD = np.divide(-asp, (1 + np.power(10, (3.9 - pH))))
        CNH2 = np.divide(1, (1 + np.power(10, (pH - 8.2))))  # ntrminal
        chrge = np.sum(np.array((CD, CE, CC, CY, CH, CK, CR, CCter, CNH2), dtype=object))
        doitisin = np.array(np.isin(np.around(chrge, 1), [0, -0])).astype(float)
        doitisin[doitisin==1] = pH
        wektor_izooelectrics = np.concatenate((wektor_izooelectrics, [doitisin]), axis=0)
        #wektor = wektor + doitisin
        pH+= 0.001
        print(pH)

    seria = []
    for seqid in range(len(wn)):
        arr = np.array(0).astype(float)
        for indx in range(len(wektor_izooelectrics)):
             arr = np.append(arr, wektor_izooelectrics[indx][seqid])
        izo_value = np.divide(np.sum(arr), len(arr[arr > 0]))
        seria.append(izo_value)
        if len(arr[arr > 0]) ==0:
            print("Zero!")
    if len(seria) == len(wn):
        print(True, len(seria), len(wn),seria)
        return seria
    else:
        print(False )


def numbering_amino_grops():
    data = pd.read_excel(r"/home/arjissuan/Desktop/MD19A/Most_active_from_MD19A.ods", engine="odf")




if __name__ == "__main__":
    #numbering_amino_grops()
    izooo()
