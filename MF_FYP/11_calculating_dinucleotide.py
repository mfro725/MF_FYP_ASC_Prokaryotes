#Script to calculate observed vs expected dinucleotide rates
#Imports
import glob 
import os
import re
from pathlib import Path
import csv


#defines paths

my_path = os.getcwd() 
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 

#list of embl files in both bacterial and archaeal folders

file_list = bacterial_genome_list + archaea_genome_list 
file_list_short = file_list[:10] #first three in list to practice running script on 


#iterates through all .embl files in both folders 

column_titles = ['accession', 'predicted_CT', 'actual_CT',	'predicted_CG', 'actual_CG', 'predicted_TG', 'actual_TG', 'predicted_CA', 'actual_CA', 'predicted_TT', 'actual_TT',	'predicted_TA', 'actual_TA','predicted_GG', 'actual_GG', 'predicted_AG', 'actual_AG', 'predicted_CC', 'actual_CC',	'predicted_GT', 'actual_GT', 'predicted_GA',	'actual_GA', 'predicted_AT', 'actual_AT',	'predicted_TC', 'actual_TC',	'predicted_AA', 'actual_AA',	'predicted_GC', 'actual_GC',	'predicted_AC', 'actual_AC']
nucleotide_content = []
csv_rows = nucleotide_content.insert(0,column_titles)

for fl in file_list:  

	genome = open(fl, "r")
	read_genome = genome.read()
	genome.close()
		
	accession = Path(fl).stem #captures accession number.embl
	print(f"accession number is {accession}")

#extracts the entire sequence

	full_seq = read_genome.rsplit(";", 1)[1] #removes annotations to give full sequence
		
#clean the sequence by removing any white spaces or numbers (using regular expression searching)

	cleaned_seq = re.sub("[^a-z]","",full_seq.lower()) #finds anything not letter, replaces with nothing in full_seq
	total_nucleotides = len(cleaned_seq)

#finding the GC content at the third position (of codon)

	nucs = "actg"
	Alldict = {}
		
		
	for nuc in nucs : 
		Alldict[nuc] = 0
			
	for nuc in nucs: 
		Alldict[nuc]= cleaned_seq.count(nuc) +Alldict[nuc]


#calculates predicted dinucleotide content from mononucleotide content
#dinucs = CT, CG, TG, CA, TT, TA, GG, AG, CC, GT, GA, AT, TC, AA, GC, AC	

#frequency of mono-nucleotides 

	A = Alldict["a"]/total_nucleotides 
	T = Alldict["t"]/total_nucleotides 
	C = Alldict["c"]/total_nucleotides 
	G = Alldict["g"]/total_nucleotides 

	
#predicted dinucleotide content calculated from mono frequencies
		
	predicted_CT = round(C*T,3)
	predicted_CG = round(C*G,3)
	predicted_TG = round(T*G,3)
	predicted_CA = round(C*A,3)
	predicted_TT = round(T*T,3)
	predicted_TA = round(T*A,3)
	predicted_GG = round(G*G,3)
	predicted_AG = round(A*G,3)
	predicted_CC = round(C*C,3)
	predicted_GT = round(G*T,3)
	predicted_GA = round(G*A,3)
	predicted_AT = round(A*T,3)
	predicted_TC = round(T*C,3)
	predicted_AA = round(A*A,3)
	predicted_GC = round(G*C,3)
	predicted_AC = round(A*C,3)

#actual dinucleotide content

#dict of dinucleotide

	dinuc_dict = {"ct" , "cg" , "tg", "ca", "tt", "ta", "gg", "ag", "cc", "gt", "ga", "at", "tc","aa", "gc", "ac"}
	Alldict3 = {}
	
	for dinuc in dinuc_dict : 
		Alldict3[dinuc] = 0
			
	for dinuc in dinuc_dict: 
		Alldict3[dinuc]= cleaned_seq.count(dinuc) +Alldict3[dinuc]
		
	total_dinucleotides = Alldict3["ct"] + Alldict3["cg"] + Alldict3["tg"] + Alldict3["ca"] + Alldict3["tt"] + Alldict3["ta"] + Alldict3["gg"] + Alldict3["ag"] + Alldict3["cc"] + Alldict3["gt"] + Alldict3["ga"] + Alldict3["at"] + Alldict3["tc"] + Alldict3["aa"] + Alldict3["gc"] + Alldict3["ac"]
	
	actual_CT = round(Alldict3["ct"]/total_dinucleotides,3)
	actual_CG = round(Alldict3["cg"]/total_dinucleotides,3)
	actual_TG = round(Alldict3["tg"]/total_dinucleotides,3)
	actual_CA = round(Alldict3["ca"]/total_dinucleotides,3)
	actual_TT = round(Alldict3["tt"]/total_dinucleotides,3)
	actual_TA = round(Alldict3["ta"]/total_dinucleotides,3)
	actual_GG = round(Alldict3["gg"]/total_dinucleotides,3)
	actual_AG = round(Alldict3["ag"]/total_dinucleotides,3)
	actual_CC = round(Alldict3["cc"]/total_dinucleotides,3)
	actual_GT = round(Alldict3["gt"]/total_dinucleotides,3)
	actual_GA = round(Alldict3["ga"]/total_dinucleotides,3)
	actual_AT = round(Alldict3["at"]/total_dinucleotides,3)
	actual_TC = round(Alldict3["tc"]/total_dinucleotides,3)
	actual_AA = round(Alldict3["aa"]/total_dinucleotides,3)
	actual_GC = round(Alldict3["gc"]/total_dinucleotides,3)
	actual_AC = round(Alldict3["ac"]/total_dinucleotides,3)
		
	
	row = [accession, predicted_CT, actual_CT, predicted_CG, actual_CG, predicted_TG, actual_TG, predicted_CA, actual_CA, predicted_TT, actual_TT,	predicted_TA, actual_TA,predicted_GG, actual_GG, predicted_AG, actual_AG, predicted_CC, actual_CC,	predicted_GT, actual_GT, predicted_GA,	actual_GA, predicted_AT, actual_AT,	predicted_TC, actual_TC,	predicted_AA, actual_AA,	predicted_GC, actual_GC,	predicted_AC, actual_AC ]

	
#puts data into list of lists, list of accessions and list of GC content percentage average per accession folder
	nucleotide_content.append(row)
	
#writes into csv

with open('11.1_dinucleotide_content.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(nucleotide_content)



