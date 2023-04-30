#script to produce table with genome master information 
import glob 
import os
import re
import fnmatch
from pathlib import Path 
import csv
from csv import writer
from csv import reader
import random


#makes list of accession number folder directories 
my_path = os.getcwd() #my path is wherever the script is sitting (Python_Learning)
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 
file_list = bacterial_genome_list + archaea_genome_list #list of embl files in both bacterial and archaeal folders
UTR_folder_path = os.path.join(my_path,'UTR_CDS')

file_list_short = file_list[:5] #first three in list to practice running script on 


#makes a list of accession folders in 3UTR file

column_titles = ['accession', 'taxa', 'translation_table', 'num_CDS','species']
accession_folder_list = []
rows_for_csv = []
rows_for_csv.insert(0,column_titles)


for fl in file_list:  #iterates through all files in both folders 
	

	row = []
#captures accession number.embl
	file_name = Path(fl).stem 
	row.append(file_name)
	print(f'Loop {file_name}')

#finds taxa of genome	
	string_directory = str(fl)
	if "archaea"in string_directory:
		 taxon = "archaea"
	else:
		taxon = "bacteria"
	row.append(taxon)

#reads EMBL file
	genome = open(fl, "r")
	read_genome = genome.read()
	genome.close()

#finds tt of genome file

	re_pattern_genome_tt_4 = 'transl_table=4'
	re_pattern_genome_tt_11 = 'transl_table=11'
	
	if_genome_tt_4 = bool(re.search(re_pattern_genome_tt_4,read_genome))
	if_genome_tt_11 = bool(re.search(re_pattern_genome_tt_11,read_genome))

	if if_genome_tt_4 == True:
		trans_table = '4'
		
	if if_genome_tt_11 == True:
		trans_table = '11'
		
	row.append(trans_table)



#appends accession folder list
	accession_folder = os.path.join(my_path,'3_UTR',file_name)
	accession_folder_list.append(accession_folder)
	
	
#counts the number of CDS 

	re_pattern_CDS = 'FT   CDS'
	num_CDS = len(re.findall(re_pattern_CDS, read_genome))
	row.append(num_CDS)
	
#finds the species of the bacteria or archea 
#\nOS   Polynucleobacter asymbioticus QLW-P1DMWA-1\n

	
	species_search = re.search('\nOS   (.*)\n', read_genome)
	species = species_search.group(1)
	row.append(species)

#add to list of lists for csv rows
	rows_for_csv.append(row)
	
	
	
# write into master table

with open('12.1_master_table.csv', 'w', newline='') as f:
	writer = csv.writer(f)
	writer.writerows(rows_for_csv)	
 
   	
