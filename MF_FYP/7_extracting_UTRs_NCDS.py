#script to capture the 3'UTR regions of non-coding sequence for Null 2 ASC rates

#imports
import glob 
import os
import re
from pathlib import Path
import csv

#defines pathways
my_path = os.getcwd() #my path is wherever the script is sitting (Python_Learning)
UTR_NCDS_path = os.path.join(my_path,'UTR_NCDS') 
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 
file_list = bacterial_genome_list + archaea_genome_list #list of embl files in both bacterial and archaeal folders

file_list_short = file_list[:3] #for practice runs

#iterates through all .embl files and makes accession folders in 3_UTR_NCDS master folder
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
	
	accession_folder = os.path.join(UTR_NCDS_path,file_name)
	
	#determines tt of genome from .embl file

	re_pattern_genome_tt_4 = 'transl_table=4'
	re_pattern_genome_tt_11 = 'transl_table=11'
	
	if_genome_tt_4 = bool(re.search(re_pattern_genome_tt_4,read_genome))
	if_genome_tt_11 = bool(re.search(re_pattern_genome_tt_11,read_genome))
	
	if if_genome_tt_4 == True:
		genome_tt = 4
	if if_genome_tt_11 == True:
		genome_tt = 11
	
	#makes accession folders 
	
	try:
		os.chdir(UTR_NCDS_path)
		os.mkdir(file_name)
		print(f"Directory {file_name} created") 
	except FileExistsError:
		print(f"Directory {file_name} already exists")


#extract the sequence

	full_seq = read_genome.rsplit(";", 1)[1] #removes annotations to give full sequence
		
#clean the sequence by removing any white spaces or numbers (using regular expression searching)

	cleaned_seq = re.sub("[^a-z]","",full_seq.lower()) #finds anything not letter, replaces with nothing in full_seq

#extract gene annotation

	just_annotations = read_genome.split("XX\nSQ" )[0]
	annotation = just_annotations.split("FT   gene")
	if len(annotation) <50 :
		annotation = just_annotations.split("FT   CDS")
			
	gene_list = annotation[1:]
		
	print(f"the number of genes is {len(annotation)}")

#for each NCDS capture 18 bases downstream of the stop codon (3UTR)

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

#capture locus tag 
		try : 
			re_pattern_2 = 'locus_tag="(.+?)"'
			all_locus_tags = re.search(re_pattern_2, gene)
			locus_tag = all_locus_tags.group(1)
		except :
			locus_tag = gene_count 
						
#capture translation table
		try: 
			re_pattern_3 = '/transl_table=(\d+)'
			all_translation_table = re.search(re_pattern_3, gene)
			translation_table = int(all_translation_table.group(1))
		except:
			translation_table = "na"
		

#checking our if it looks like non-coding sequence with two defined assumptions:
#sequence that is not a multiple of three 
#sequence that does not end with a primary stop codon

	#multiple of 3 long? 
		
		gene_length = len(gene_sequence)
		
		if gene_length%3 != 0: 
			gene_length_m3 = "not_multiple_3"
		else:
			gene_length_m3 = "multiple_3"
		
    		
	#end with a stop?
		
		suffixes = ('tag','taa','tga')
		has_stop_codon = bool(gene_sequence.endswith(suffixes))
		
		if has_stop_codon == True:
			stop_codon = "seq_has_sc"
		else:
			stop_codon = "seq_has_no_sc"

	#does it not have any in frame stops?
		 
		re_pattern_4 = "tag" or "taa" or "tga"
		gene_sequence_no_SC = gene_sequence[:-3]
		IF_stops = bool(re.search(re_pattern_4,gene_sequence_no_SC))
		if IF_stops == "true":
			in_frame_stops = "yes"
		else:
			in_frame_stops = "no"

#write sequence into file named after locus tag to correct folder

		if gene_length_m3 == "not_multiple_3" and stop_codon == "seq_has_no_sc":
			non_coding = 'yes'
			UTR_file_name = str(locus_tag) + "_non_coding"+ "_3UTR" + ".fta"
			
			UTR_NCDS_file_path = os.path.join(UTR_NCDS_path, file_name, UTR_file_name)
			with open(UTR_NCDS_file_path , "w") as f:
				f.write(f">{file_name};{locus_tag};tt = {translation_table};tt_genome = {genome_tt}; comp={comp};len={len(three_UTR)}; taxon = {taxon}\n{three_UTR}")
				f.close() 
		
		else:
			non_coding = 'no'
			continue 


