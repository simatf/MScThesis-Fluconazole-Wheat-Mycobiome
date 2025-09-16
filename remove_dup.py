import pandas as pd

#--------------------------------------------

df = pd.read_csv('../list_all_mod_comb.csv', usecols= [1, 2, 3, 4, 10, 14]) # 5, 6 # 7, 8,

df['Putative species'] = df['Penicillium grouping'].where(df['Penicillium grouping'].notna(), df['Putative species'])

df = df[~df['Putative species'].str.lower().str.startswith("unknown")]

df.drop(columns=['Penicillium grouping'], inplace= True)

df.drop_duplicates(inplace= True)

df.drop(columns=['Dilution'], inplace= True) # , 'Date', 'Days post initial isolation' #  'Morphotype', 'Genotype'

df.to_csv('shannon_rem_dup_no_isolates.csv', index= False)