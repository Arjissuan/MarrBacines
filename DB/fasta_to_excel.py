from Bio import SeqIO
import pandas as pd
import os
def converter(file_name):
    record = SeqIO.parse(file_name, 'fasta')
    row = []
    for item in record:
        # print(item.__dict__)
        row.append([item.seq, item.id, item.name, item.description])
    df = pd.DataFrame(data=row, columns=['sequence', 'id', 'name', 'description'])
    seq = df['sequence']
    for i in range(len(seq)):
        seq.iloc[i] = ''.join(seq.iloc[i])
    df['sequence'] = seq

    id = df['id']
    name = df['name']
    description = df['description']
    for i in range(len(seq)):
        # description.iloc[i] = description.iloc[i].replace(id.iloc[i], '')
        id_name = id.iloc[i].split('|')
        id.iloc[i] = id_name[0]
        name.iloc[i] = id_name[1]
    df['id'] = id
    df['name'] = name
    df['description'] = description

    return df

if __name__ == '__main__':
    file_name = input('Provide name of the fasta file with extension: ')
    name = input('Provide name of the out file without extension: ')
    df = converter(file_name)
    df.to_excel(os.path.join(os.getcwd(),f'{name}.xlsx'))