#script to produce expected ASC rates from OOF sequence (Null 3) 
#imports
import glob 
import os
import re
import fnmatch
from pathlib import Path 
import openpyxl
import csv
import pandas as pd


#makes list of accession number folder directories 
my_path = os.getcwd() #my path is wherever the script is sitting (Python_Learning)
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 
file_list = bacterial_genome_list + archaea_genome_list #list of embl files in both bacterial and archaeal folders
UTR_folder_path = os.path.join(my_path,'UTR_OOF')

file_list_short = file_list[4:5] #for practice runs


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
	
	accession_folder = os.path.join(my_path,'UTR_OOF',file_name)
	accession_folder_list.append(accession_folder)
	
accession_folder_list_short = accession_folder_list[:5] #short list for test runs


#make excel file

data_columns = ['accession', 'tot_UTRs_OOF', 'tot_seq_OOF', 'num_of_tag_OOF', 'num_of_tga_OOF', 'num_of_taa_OOF', 'total_number_of_ASC_OOF', 'frequency_tag_OOF', 'frequency_tga_OOF', 'frequency_taa_OOF' , 'frequency_ASC_OOF']

data_row_lists = []

csv_rows = data_row_lists.insert(0,data_columns)

#loops through accession folders to count stops 

for flder in accession_folder_list: 
	
	UTR_files = glob.glob(os.path.join(flder,'*+2_shifted_UTR.fta')) #list of 3UTR files in accession folder
	
	accession_number = os.path.basename(flder) #accession number folder 
	
	print (f'accession name is {accession_number}')
	
	total_tag_c = 0 #total tag in protein coding coding (per accession)
	total_tga_c  = 0 
	total_taa_c = 0 

	
	for fl in UTR_files:
	
		tot_seq = 18*len(UTR_files)
		tot_UTR = len(UTR_files)
	
		UTR = open(fl, "r")
		read_UTR = UTR.read()
		UTR.close()
		
#take 18 bases from 3UTR file and repackage into codons

		UTR_seq = read_UTR.rsplit("\n", 1)[1]
		cdns =[]

		for i in range(0, len(UTR_seq), 3) :
			cdns.append(UTR_seq[i:i + 3])
	
	
#check for stops
		
		r_tag = re.compile('tag')
		r_tga = re.compile('tga')
		r_taa = re.compile('taa')
		
		re_pattern_tt_11 = 'tt = 11'
		re_pattern_tt_4 = 'tt = 4'
		re_pattern_tt_na = 'tt = na'		
		
		if_tt_11 = bool(re.search(re_pattern_tt_11,read_UTR))
		if_tt_4 = bool(re.search(re_pattern_tt_4,read_UTR))
		if_tt_na = bool(re.search(re_pattern_tt_na,read_UTR))
		

		#tt=11 or tt=4 protein coding totals
		
		if if_tt_11 == True :
			
			numstops_tag = len(list(filter(r_tag.match, cdns)))
			total_tag_c = total_tag_c + numstops_tag
			
			numstops_tga = len(list(filter(r_tga.match, cdns)))
			total_tga_c = total_tga_c + numstops_tga

			numstops_taa = len(list(filter(r_taa.match, cdns)))
			total_taa_c = total_taa_c + numstops_taa

		if if_tt_4 == True : #tga not recognised as a stop in tt = 4
		
			numstops_tag = len(list(filter(r_tag.match, cdns)))
			total_tag_c = total_tag_c + numstops_tag
		
			numstops_taa = len(list(filter(r_taa.match, cdns)))
			total_taa_c = total_taa_c + numstops_taa
			
	numstops_c = total_tag_c + total_tga_c + total_taa_c 
	
	#frequency of of sc per codon 
	frequency_tag = round((total_tag_c /(tot_UTR*6)),4)
	frequency_tga = round((total_tga_c / (tot_UTR*6)),4)
	frequency_taa = round((total_taa_c / (tot_UTR*6)),4)
	frequency_all = round((numstops_c/ (tot_UTR*6)),4)
	
	#makes a list of data to go into row in excel 
	accession_data_list_row = [accession_number,tot_UTR, tot_seq, total_tag_c,total_tga_c,total_taa_c, numstops_c,frequency_tag, frequency_tga, frequency_taa, frequency_all] #lemgth of number of
	#appends list of lists (lists of data to go into rows)
	data_row_lists.append(accession_data_list_row)


#proportion to the number of genes/size of genomes

#write list of lists into csv

with open('10.1_null3_ASC_rates.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(data_row_lists)

	

