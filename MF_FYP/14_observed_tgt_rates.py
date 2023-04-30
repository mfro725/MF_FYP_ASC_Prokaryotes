#script to produce expected ASC rates from 3'UTR of NCDS (Null 2)

import glob 
import os
import re
import fnmatch
from pathlib import Path 
import openpyxl
import csv


#makes list of accession number folder directories 
my_path = os.getcwd() #my path is wherever the script is sitting (Python_Learning)
UTR_NCDS_path = os.path.join(my_path,'UTR_CDS') 
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 
file_list = bacterial_genome_list + archaea_genome_list #list of embl files in both bacterial and archaeal folders


file_list_short = file_list[:6] #first three in list to practice running script on 


#makes a list of accession folders in 3UTR file


accession_folder_list = []

for fl in file_list:  #iterates through all files in both folders 
	if fl.find("archaea"):
		 taxon = "archaea"
	else:
		taxon = "bacteria"
	genome = open(fl, "r")
	read_genome = genome.read()
	genome.close()
		
	file_name = Path(fl).stem #captures accession number.embl
	#print(f"accession number is {file_name}")
	
	accession_folder = os.path.join(my_path,'UTR_CDS',file_name)
	accession_folder_list.append(accession_folder)
	
accession_folder_list_short = accession_folder_list[:5] #short list for test runs


#make excel file

data_columns = ['accession', 'tot_UTRs_coding', 'tot_seq_coding', 'num_of_tgt_coding', 'frequency_tgt_coding']

data_row_lists = []

csv_rows = data_row_lists.insert(0,data_columns)



#loops through accession folders to count stops 

for flder in accession_folder_list: 

	
	UTR_files = glob.glob(os.path.join(flder,'*_3UTR.fta')) #list of 3UTR files in accession folder
	
	accession_number = os.path.basename(flder) #accession number folder 
	print (f'accession is {accession_number}')
	total_tgt_c = 0 #total tgt in protein coding coding (per accession)
	

#frequency of of stops per codon (number of stops / total number of UTR regions * 6)

	tot_UTR = len(UTR_files)
	tot_seq = 18*len(UTR_files)
	
	for fl in UTR_files:

	
		UTR = open(fl, "r")
		read_UTR = UTR.read()
		UTR.close()
		
#take 18 bases from 3UTR file and repackage into codons

		UTR_seq = read_UTR.rsplit("\n", 1)[1]
		cdns =[]

		for i in range(0, len(UTR_seq), 3) :
			cdns.append(UTR_seq[i:i + 3])
	
	
#check for stops
		
		r_tgt = re.compile('tgt')
		
		re_pattern_tt_11 = 'tt = 11'
		re_pattern_tt_4 = 'tt = 4'
		re_pattern_tt_na = 'tt = na'		
		
		if_tt_11 = bool(re.search(re_pattern_tt_11,read_UTR))
		if_tt_4 = bool(re.search(re_pattern_tt_4,read_UTR))
		if_tt_na = bool(re.search(re_pattern_tt_na,read_UTR))
		
		#tt=na non protein coding totals 
		if if_tt_na == True :
		
			numstops_tgt = len(list(filter(r_tgt.match, cdns)))
			total_tgt_c = total_tgt_c + numstops_tgt
			
		
		#tt=11 or tt=4 protein coding totals
		
		if if_tt_11 == True :
			
			numstops_tgt = len(list(filter(r_tgt.match, cdns)))
			total_tgt_c = total_tgt_c + numstops_tgt
			

	
	#frequency of of sc per codon 
	if tot_UTR == 0:
		frequency_tgt = 'NA'
		
	else: 
		frequency_tgt = round((total_tgt_c /(tot_UTR*6)),4)

		
	#makes a list of data to go into row in excel 
	accession_data_list_row = [accession_number,tot_UTR, tot_seq, total_tgt_c, frequency_tgt] #lemgth of number of
	
	#appends list of lists (lists of data to go into rows)
	data_row_lists.append(accession_data_list_row)


#proportion to the number of genes/size of genomes

#write list of lists into csv

with open('tgt_observed_test.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(data_row_lists)

