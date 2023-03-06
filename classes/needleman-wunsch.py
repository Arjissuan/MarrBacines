import numpy as np
import pandas as pd



def nw1(seq1, seq2, match=3, missmatch=-1, gap=-3):
    prep = np.zeros([len(seq2)+1, len(seq1)+1])
    prep[0, 1:] = np.linspace(gap-1, gap-len(seq1), len(seq1))
    prep[1:, 0] = np.linspace(gap-1, gap-len(seq2), len(seq2))
    DT = pd.DataFrame(
        prep,
        columns=list(nucl for nucl in ("s" + seq1)),
        index=list(nucle for nucle in ("s" + seq2))
    )

    DTwun = DT.copy()
    DTpath = DT.copy()

    for col in DT.columns:
        if col != "s":
            for indx in DT.index:
                if indx != "s":
                    if col == indx:
                        DT.at[indx, col] = match
                    else:
                        DT.at[indx, col] = missmatch

    #print(DTwun, "DTwun\n")
    print(DT, "DT (these are small s values)\n")
    #needlmen_wunsch algorithm, making socre for each match
    for col in range(1, len(DT.columns)):
        for indx in range(1, len(DT.index)):
            a = int(DTwun.iat[indx - 1, col - 1]) #big S values
            b = int(DT.iat[indx, col]) #small s value
            c = int(DTwun.iat[indx - 1, col]) #up
            d = int(DTwun.iat[indx, col-1]) #left
            DTwun.iat[indx, col] = np.max(  #makes scores
                (
                    np.sum((a, b)),
                    np.sum((d, -4)),
                    np.sum((c, -4))
                )

            )
            #makes path
            if np.sum((a, b)) > np.sum((d, -4)) and np.sum((a, b)) > np.sum((c, -4)):
                DTpath.iat[indx, col] = 0
            elif np.sum((c, -4)) > np.sum((a,b)) and np.sum((c,-4)) > np.sum((d,-4)):
                DTpath.iat[indx, col] = -b
            elif np.sum((d, -4)) > np.sum((a, b)) and np.sum((d, -4)) > np.sum((c, -4)):
                DTpath.iat[indx, col] = b

    print(DTwun, "\DTwun after (these are big S values)")
    print(DTpath, "\DTpath (these are small t values)")
    Mwun = DTwun.to_numpy()
    Mpath = DTpath.to_numpy()
    i = np.argmax(Mwun[:,-1:])
    j = len(Mpath[0,:])-1
    nx = []
    ny = []

    #using previosly made path
    #this needs to be revorked!!!
    while j > 0 or i > 0:
        #print(i,j)
        if Mpath[i,j] == 0:
            ny.append(DT.index[i])
            nx.append(DT.columns[j])
            i+=-1
        elif Mpath[i,j] > 0:
            ny.append(DT.index[i])
            nx.append(DT.columns[j])
            nx.reverse()
            nx.append("-")
            nx.reverse()
            i+=-1
            j+=1
        elif Mpath[i,j] < 0:
            ny.append(DT.index[i])
            nx.append(DT.columns[j])
            nx.reverse()
            nx.append("-")
            nx.reverse()

        j+=-1
    if "-" in nx:
        ny.append("-")
    #print(nx, ny)
    nu = []
    for index, value in enumerate(seq2[::-1]):
        #print(index, value)
        ny[index] = value
        if nx[index] == "-":
            nu.append(' ')
        elif ny[index] == nx[index]:
            nu.append('|')
        elif ny[index] != nx[index]:
            nu.append('.')
        if "-" == nx[0] and index == len(seq2)-1:
            nu.append(" ")


    return "{}\n{}\n{}\n".format(" ".join(nx[::-1]), " ".join(nu[::-1]), " ".join(ny[::-1]))



if __name__ == "__main__":
    s1 = "GGACTG"
    s2 = "GGATGC"
    se1 = 'AGACAA'
    se2 = 'CGAACT'
    print(nw1(seq1=se1, seq2=se2))

#do poprawy cale przedstawienie aligmentu, dziala na zlym motywie. Ramka odczytu przesuwa sie wraz z miejsce rozpoczecia na macierzy, a nie jak jest obecnie.