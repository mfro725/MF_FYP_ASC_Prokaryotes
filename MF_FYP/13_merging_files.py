#merges TAG, TGA and TAA simulation rates and GC3 content files
import pandas as pd
import csv
  
# reading two csv files
data1 = pd.read_csv('3.1_null1_tag_rates_1.csv')
data2 = pd.read_csv('4.1_null1_tga_rates_1.csv')


# using merge function by setting how='inner'
output1 = pd.merge(data1, data2,
                   on='accession', 
                   how='inner')

# writing result into csv file

output1.to_csv('13.1.1_merged_simulation_rates.csv', index=False)
print(f'merger 1 done')

data3 = pd.read_csv('13.1.1_merged_simulation_rates.csv')
data4 = pd.read_csv('5.1_null1_taa_rates_1.csv')

output2 = pd.merge(data3, data4,
                   on='accession', 
                   how='inner')

output2.to_csv('13.1.2_merged_simulation_rates.csv', index=False)
print(f'merger 2 done')

data5 = pd.read_csv('13.1.2_merged_simulation_rates.csv')
data6 = pd.read_csv('6.1_GC3_content_per_genome.csv')

output3 = pd.merge(data5, data6,
                   on='accession', 
                   how='inner')

output3.to_csv('13.1_merged_simulation_rates.csv', index=False)
print(f'merger of simulation files complete')



#merges observed ASC rates file with rates produced from NCDS (null 2)
  
# reading two csv files
data7 = pd.read_csv('2.1_observed_ASC_rates.csv')
data8 = pd.read_csv('8.1_null2_ASC_rates.csv')


  
# using merge function by setting how='inner'
output4 = pd.merge(data7, data8,
                   on='accession', 
                   how='inner')

# writing result into csv file

output4.to_csv('13.2_merged_NCDS_CDS_rates.csv', index=False)








'''

#merges observed ASC rates with rates produced from OOF seq (null 3)

# reading two csv files
data9 = pd.read_csv('2.1_observed_ASC_rates.csv')
data10 = pd.read_csv('10_null3_ASC_rates.csv')


  
# using merge function by setting how='inner'
output5 = pd.merge(data9, data10,
                   on='accession', 
                   how='inner')

# writing result into csv file

output5.to_csv('13.3_merged_OOF_CDS_rates.csv', index=False)

'''





