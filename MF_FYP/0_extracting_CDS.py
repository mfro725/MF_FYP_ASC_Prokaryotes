#Script to extract protein coding sequence from folder of genome .embl files
import glob 
import os
import re
from pathlib import Path
#can make folders at level of directory 
#defines genome file pathway of both bacterial and archaea folders and creates list of embl files
my_path = os.getcwd() #my path is wherever the script is sitting 
bacterial_folder_path = os.path.join(my_path,'bacterial_filtered_genomes') 
archaea_folder_path = os.path.join(my_path,'archaea_filtered_genomes')
bacterial_embl_path = os.path.join(bacterial_folder_path,"*.embl")
archaea_embl_path= os.path.join(archaea_folder_path,"*.embl")
bacterial_genome_list = glob.glob(bacterial_embl_path) 
archaea_genome_list = glob.glob(archaea_embl_path) 
CDS_folder_path = os.path.join(my_path,'CDS') 
file_list = bacterial_genome_list + archaea_genome_list #list of embl files in both bacterial and archaeal folders

file_list_short = file_list[:3] #first three in list to practice running script on 

for fl in file_list:  #iterates through all files in both folders 
	genome = open(fl, "r")
	read_genome = genome.read()
	genome.close()
		
	file_name = Path(fl).stem #captures accession number.embl
	print(f"accession number is {file_name}")
	
	
	accession_folder = os.path.join(my_path,'CDS',file_name)
	
	try:
		os.mkdir(accession_folder)
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


#for each gene capture gene position, extract gene sequence from full sequence

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
			translation_table = all_translation_table.group(1)
		except:
			translation_table = "na"
		
	
#checking does this look like a 'good sequence':

	#Is it a multiple of 3 long? 
		
		gene_length = len(gene_sequence)
		
		try:
			gene_length % 3 == 0
			gene_length_m3 = "multiple of 3"
		except:
			continue 
    		
	#Does it end with a stop (TAG, TAA and TGA)
	
		try: 
			
			suffixes = ('tag','taa','tga')
			has_stop_codon = gene_sequence.endswith(suffixes) #was has_stop_codon = gene_sequence.endswith('tag'or'taa'or'tga')
			stop_codon = gene_sequence[-3:]
		except:
			stop_codon = "na"

	#Does it not have any in frame stops? If yes, then read through.
		 
		re_pattern_4 = "tag" or "taa" or "tga"
		gene_sequence_no_SC = gene_sequence[:-3]
		IF_stops = bool(re.search(re_pattern_4,gene_sequence_no_SC))
		if IF_stops == "true":
			in_frame_stops = "yes"
		else:
			in_frame_stops = "no"

#write sequence into file named after locus tag to correct folder
		gene_file_name = str(locus_tag) + ".fta"
		file_path = os.path.join(accession_folder,gene_file_name)
		with open(file_path, "w") as f:
			f.write(f">{locus_tag};tt = {translation_table};SC = {stop_codon};IF_SC = {in_frame_stops};comp={comp};len={len(gene_sequence)}\n{gene_sequence}")
			f.close()
				


