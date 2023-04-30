#merges TGT and TGG observed and expected rates from null 1 and 2 for R analysis 

import pandas as pd
import csv

df1 = pd.read_csv('19.1_observed_tgt_rates.csv')
df2 = pd.read_csv('20.1_observed_tgg_rates.csv')
df3 = pd.read_csv('21.1_NCDS_tgt_rates.csv')
df4 = pd.read_csv('22.1_NCDS_tgg_rates.csv')
df5 = pd.read_csv('23.1_OOF_tgt_rates.csv')
df6 = pd.read_csv('24.1_OOF_tgg_rates.csv')

#files = ['19.1_observed_tgt_rates.csv', '20.1_observed_tgg_rates.csv', '21.1_NCDS_tgt_rates.csv']


df_main = df1.merge(df2,on='accession', how = 'inner').merge(df3,on='accession', how = 'inner').merge(df4,on='accession', how = 'inner').merge(df5,on='accession', how = 'inner').merge(df6,on='accession', how = 'inner')

df_main.to_csv('25.1_tgt_tgg_data.csv', index = False)