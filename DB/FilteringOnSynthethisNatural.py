import pandas as pd
df = pd.read_excel('general_amps.xlsx')
q = df['Source'].str.lower().str.contains('synthetic', 'synthesized')
print(len(q))
indx = q[q == True].index
df.loc[indx, :].to_excel('Synthethised.xlsx')
indx_ = q[q == False].index
df.loc[indx_, :].to_excel('Natural.xlsx')