import numpy as np
import pandas as pd
def data_cleaning(csv_file):
    with open(csv_file, 'r') as my_file:
        for line in my_file.readlines():
            tabnum = np.sum(list(map(lambda x: (1 if x == '\t' else 0), line)))
            if tabnum != 12:
                print(line)

def great_replacment(dataframe, columns):
    for col in columns:
        purged = dataframe.loc[:, col].str.replace('[', '', ).str.replace(']', '').str.replace('\'', '')
        dataframe.loc[:, col] = purged
    return dataframe

def grand_merge(dataframe):
    df = dataframe
    merger = lambda x: (f'({df.iloc[x, 6]})({df.iloc[x, 7]})({df.iloc[x, 8]})({df.iloc[x, 9]})'
                        f'({df.iloc[x, 10]})({df.iloc[x, 10]})({df.iloc[x, 11]})({df.iloc[x, 12]})')
    merge_col = list(map(merger, range(len(dataframe.index))))
    df['origin'] = merge_col
    return df



if __name__ == '__main__':
    # df = pd.read_csv('/home/arjissuan/Desktop/DBAASP.csv', delimiter='\t')
    # df.to_excel('/home/arjissuan/Desktop/DBAASP.xlsx')
    # data_cleaning('/home/arjissuan/Desktop/DBAASP.csv')
    df = pd.read_excel('/home/arjissuan/Desktop/AMP-databases/DBAASP.xlsx', index_col=0)
    columns_to_purge = ['targeObjects', 'targetGroups', 'kingdom', 'subkingdom', 'source', 'plasmid', 'gene']
    upd_df = great_replacment(df, columns_to_purge)
    new_df = grand_merge(upd_df)
    new_df.to_excel('/home/arjissuan/Desktop/AMP-databases/DBAASP.xlsx',)
