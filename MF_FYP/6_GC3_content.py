#script to calculate GC content of genomes. GC3% per accession number.
import glob 
import os
import re
from pathlib import Path
import numpy as np
import csv

#create list of all .fta gene files 

my_path = os.getcwd() #my path is wherever the script is sitting (Python_Learning)
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_all_path = os.path.join(bacterial_folder_path,"*")
archaea_all_path= os.path.join(archaea_folder_path,"*")
bacterial_all_list = glob.glob(bacterial_all_path) 
archaea_all_list = glob.glob(archaea_all_path) 
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path)
CDS_folder_path = os.path.join(my_path,'CDS') 
CDS_all_list = os.path.join(CDS_folder_path,"*") 

 #loop to make list of all accession folders 
 
genome_folder_list = glob.glob(CDS_all_list)

genome_folder_list_short = genome_folder_list[610:620] #for test runs


#for accession folder in genome_folder_list makes list of just gene files (not _3UTR)

column_titles = ['accession','GC3']

GC_percentage = []

csv_rows = GC_percentage.insert(0,column_titles)

for flder in genome_folder_list:
	
	total_gc_content = 0
	accession_number = os.path.basename(flder)
	

	fasta_path = os.path.join(flder,"*.fta")
	fasta_files = glob.glob(fasta_path)
	
	CDS_list = []
	for element in fasta_files:
		CDS_list.append(element)
	
	
	for fl in CDS_list: 
		
		locus_tag = os.path.basename(fl)
	
		gene = open(fl, "r")
		read_gene = gene.read() #read faster file
		gene.close()
		
		#extracts gene sequence from .fta file and packages into codons 
		gene_sequence = read_gene.rsplit("\n", 1)[1]
	
		
	#finding the GC content at the third position (of codon)

		nucs = "actg"
		s = gene_sequence[::3]
		Alldict = {}
		
		
		for nuc in nucs : 
			Alldict[nuc] = 0
			
		for nuc in nucs: 
			Alldict[nuc]= s.count(nuc) +Alldict[nuc]

		CountAllGC = Alldict["g"] + Alldict["c"]
		TotalNucleotides =  Alldict["g"] + Alldict["c"] +  Alldict["a"] + Alldict["t"]
		PercentageGCContent = (CountAllGC / TotalNucleotides)*100
		
		
		total_gc_content = total_gc_content + PercentageGCContent
		
	total_gc_content_pergentage_accession = round(total_gc_content /len(CDS_list),2)
	
	print(f'The total GC content for accession {accession_number} is {total_gc_content_pergentage_accession}%')

#puts data into list of lists, list of accessions and list of GC content percentage average per accession folder
	data_list = [accession_number,total_gc_content_pergentage_accession]
	GC_percentage.append(data_list)

#export to excel
	
	with open ("6.1_GC3_content_per_genome.csv",'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerows(GC_percentage)












