#script to produce distribution of tga rates for 1000 simulated 3'UTR sequence

#imports
import glob 
import os
import re
import fnmatch
from pathlib import Path 
import csv
from csv import writer
from csv import reader
import numpy as np
import random

#defines pathways

my_path = os.getcwd() 
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 
file_list = bacterial_genome_list + archaea_genome_list #list of embl files 
UTR_folder_path = os.path.join(my_path,'UTR_CDS')

file_list_short = file_list[:3] #short list for test runs

##makes a list of accession folders in 3UTR file

accession_folder_list = []

accession_tt = {}

for fl in file_list:  
		
	genome = open(fl, "r")
	read_genome = genome.read()
	genome.close()

#determines tt of genome from .embl file
	re_pattern_genome_tt_4 = 'transl_table=4'
	re_pattern_genome_tt_11 = 'transl_table=11'
	
	if_genome_tt_4 = bool(re.search(re_pattern_genome_tt_4,read_genome))
	if_genome_tt_11 = bool(re.search(re_pattern_genome_tt_11,read_genome))
	
	file_name = Path(fl).stem 
	
	if if_genome_tt_4 == True:
		accession_tt[file_name] = 4
	if if_genome_tt_11 == True:
		accession_tt[file_name] = 11


#captures accession number and appends accession list with folde directories
		
	
	accession_folder = os.path.join(my_path,'UTR_CDS',file_name)
	accession_folder_list.append(accession_folder)
	
accession_folder_list_short = accession_folder_list[:4] #short list for test runs

cnt = 0 #count for tracking running

#master list of lists to be written into csv file


number_tga_c_randomised_list = [] #list of rows (1000 columns)


column_name_list = [] #list of titles for columns
	
for randomisation_number in range(1,1001):
	column_name = str(randomisation_number) + '_shuffle_c'
	column_name_list.append(column_name)
		


column_name_list.insert(0,'Z_score')
column_name_list.insert(0,'mean')
column_name_list.insert(0,'st_dev')
column_name_list.insert(0,'observed_tga')
column_name_list.insert(0,'tot_coding')
column_name_list.insert(0,'accession')


rows = []

rows.insert(0,column_name_list)

#reads in actual stop codon rates 2.1_observed_ASC_rates.csv file to use to calc Z value. 
#column 5 = tga

list_actual_tga_rates = []

with open('2.1_observed_ASC_rates.csv') as fl_obj:
	reader_obj = csv.reader(fl_obj)
	for row in reader_obj:
		list_actual_tga_rates.append(row[4])
		
del list_actual_tga_rates [0] #removes column title


#loops through files in accession folder to total tga in 3'UTR regions for each randomisation
#calculates mean, standard deviation and Z scores using observed rates

for flder in accession_folder_list:

	UTR_files = glob.glob(os.path.join(flder,'*_3UTR.fta')) 

    #for tracking run
	accession_number = os.path.basename(flder) 
	cnt = cnt + 1
	print(f'accession is {accession_number} : {cnt} of {len(accession_folder_list)}')
	
	accession_row_list_c = []
	


	#shuffled 1000 times
	
	for shuffle in range(1000): 
		
		if shuffle%5 == 0 :
			print(f'shuffle {shuffle}')
		
		total_tga_all_files_c = 0
		
		
		for fl_1 in UTR_files:
			
			number_tga_c_shuffle_i = 0
	
	#opens, reads, closes UTR file in accession folder 
			UTR = open(fl_1, "r")
			read_UTR = UTR.read()
			UTR.close()
		
	#extracts UTR sequence 
			UTR_seq = read_UTR.rsplit("\n", 1)[1]
		
	#shuffles UTR sequence once 
			UTR_seq_list = list(UTR_seq)
			random.shuffle(UTR_seq_list)
			UTR_shuffled = ''.join(UTR_seq_list) 
			
	#package UTRs into codons
			cdns =[]
			for i in range(0, len(UTR_shuffled), 3) :
				cdns.append(UTR_shuffled[i:i + 3])
			
	#count tga 
			r_tag = re.compile('tag')
			r_tga = re.compile('tga')
			r_taa = re.compile('taa')
		
			re_pattern_tt_11 = 'tt = 11'
			re_pattern_tt_4 = 'tt = 4'	
		
			if_tt_11 = bool(re.search(re_pattern_tt_11,read_UTR))
			if_tt_4 = bool(re.search(re_pattern_tt_4,read_UTR))
		
				
	#tt_11 count tga only. tt_4 does not recognise tga as a stop.
			
			
			if if_tt_11 == True :
			
				numstops_tga = len(list(filter(r_tga.match, cdns)))
				number_tga_c_shuffle_i = number_tga_c_shuffle_i + numstops_tga
		

			total_tga_all_files_c = total_tga_all_files_c + number_tga_c_shuffle_i
		
			tot_coding = len(UTR_files)
		
		accession_row_list_c.append(total_tga_all_files_c)
		
		
	#calculate mean and standard deviation of null distribution
	
	st_dev = round((np.std(accession_row_list_c)),2)
	shuffle_mean = round((np.mean(accession_row_list_c)),2)
	
	#calculate Z score 
	
	column_value = cnt - 1
	observed_tga = list_actual_tga_rates[column_value]
	top = float(observed_tga) - (shuffle_mean)
	Z_score = round((top / st_dev),2)

	#calculates total UTR regions in accession folder for frequency calc
	
	
	
	rows_for_csv = [accession_number, tot_coding, observed_tga, st_dev, shuffle_mean, Z_score, *accession_row_list_c]
		
	rows.append(rows_for_csv) #appends list of lists for csv file

# write into file 

with open('4.1_null1_tga_rates_1.csv', 'w', newline='') as f:
   	writer = csv.writer(f)
   	writer.writerows(rows)
