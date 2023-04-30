#script to capture the 18 bases (6 codons) downstream of the final stop codon (3'UTRs)
#puts 3UTRs into files named : locus_tag_3UTR in appropriate accession number folders
#puts accession folders in master 'UTRs' folder

#imports
import glob 
import os
import re
from pathlib import Path
import csv


#defines genome file pathways and creates list of embl files

my_path = os.getcwd() 
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 

file_list = bacterial_genome_list + archaea_genome_list 

file_list_short = file_list[:3] #for test runs

#iterates through all files in bacterial_filtered_genomes and archea_filtered_genomes

for fl in file_list:  

	if fl.find("archaea"):
		 taxon = "archaea"
	else:
		taxon = "bacteria"
	genome = open(fl, "r")
	read_genome = genome.read()
	genome.close()
		
	file_name = Path(fl).stem #captures accession number.embl
	print(f"accession number is {file_name}")
	
	accession_folder = os.path.join(my_path,'UTR_CDS',file_name)
	
	try:
		os.mkdir(accession_folder)
		print(f"Directory {accession_folder} created") 
	
	except FileExistsError:
		print(f"Directory {file_name} already exists")

		
#extracts the CDS 

	full_seq = read_genome.rsplit(";", 1)[1] #removes annotations to give full sequence
		
#cleans the sequence by removing any white spaces or numbers (using regular expression searching)

	cleaned_seq = re.sub("[^a-z]","",full_seq.lower()) #finds anything not letter, replaces with nothing in full_seq

#extracts gene annotation

	just_annotations = read_genome.split("XX\nSQ" )[0]
	annotation = just_annotations.split("FT   gene")
	if len(annotation) <50 :
		annotation = just_annotations.split("FT   CDS")
			
	gene_list = annotation[1:]
		
	print(f"the number of genes is {len(annotation)}")

#captures 18 bases downstream of the primary stop (3UTR) for each gene

	complementary_bases = {"a":"t", "t":"a", "c":"g", "g":"c"}
	gene_count = 0 #if you can't find locus tag it will add to count 
	
	for gene in gene_list : 
		gene_count = gene_count + 1
		gene_bits = gene.split("\n")
		top_line = gene_bits[0].strip()
		try:
			re_pattern = "(\d+)\.\.(\d+)"  #() captures group, group (1) is the first pattern (diget), group (2) is the second and group(0) is complete capture
			coordinates_both = re.search(re_pattern,top_line)
			unshifted_start = int(coordinates_both.group(1))
			shifted_start = unshifted_start - 1  #shifts because EMBL file starts at 1 and python begins at 0
			finish = int(coordinates_both.group(2))
			gene_sequence = (cleaned_seq[shifted_start:finish])
			end_3UTR = finish + 18
			three_UTR = (cleaned_seq[finish:end_3UTR])
			
		except:
			continue
			
#if complement, reverse sequence and filp bases

		comp = "no"  #variable for each gene to state if complement or not 
		if "complement" in top_line:
			comp = "yes" 
			gene_sequence = gene_sequence[::-1]
			reversed_seq = ""
			for base in gene_sequence:
				if base in complementary_bases: 
					reversed_seq = reversed_seq + complementary_bases[base] 
				else :
					reversed_seq = reversed_seq + base 
			gene_sequence = reversed_seq 
			

#captures locus tag 

		try : 
			re_pattern_2 = 'locus_tag="(.+?)"'
			all_locus_tags = re.search(re_pattern_2, gene)
			locus_tag = all_locus_tags.group(1)
		except :
			locus_tag = gene_count 
						
#captures translation table

		try: 
			re_pattern_3 = '/transl_table=(\d+)'
			all_translation_table = re.search(re_pattern_3, gene)
			translation_table = int(all_translation_table.group(1))
		except:
			translation_table = "na"
		
	
#checking does this look like a protein-coding:

	#Is it a multiple of 3 long? 
	
		protein_coding = 1
		gene_length = len(gene_sequence)
		
		try:
			gene_length % 3 == 0
			gene_length_m3 = "multiple of 3"
		except:
			protein_coding = 0
    		
	#Does it end with a stop (TAG, TAA and TGA)

		try: 
			has_stop_codon = gene_sequence.endswith("tag"or"taa"or"tga")
			stop_codon = gene_sequence[-3:]
		except:
			protein_coding = 0 

	#Does it not have any in frame stops? 
		 
		re_pattern_4 = "tag" or "taa" 
		re_pattern_11 = "tag" or "taa" or "tga"
		gene_sequence_no_SC = gene_sequence[:-3]
		if translation_table == 4:
			IF_stops = bool(re.search(re_pattern_4,gene_sequence_no_SC))
		elif translation_table == 11: 
			IF_stops = bool(re.search(re_pattern_11,gene_sequence_no_SC))
		else: 
			protein_coding = 0 
		if IF_stops == "true":
			protein_coding = 0
		

#writes into files and puts into accession folder in UTRs master folder
		UTR_file_name = str(locus_tag) + "UTR_CDS" + ".fta"
		if protein_coding == 0 : 
			CDS = "no"
		else : 
			CDS = "yes"
			
		UTR_file_path = os.path.join(accession_folder,UTR_file_name)
		with open(UTR_file_path, "w") as f:
			f.write(f">{file_name};{locus_tag};tt = {translation_table};SC = {stop_codon};Prot = {CDS};comp={comp};len={len(three_UTR)}; taxon = {taxon}\n{three_UTR}")
			f.close()
				



