# MF_FYP_ASC_Prokaryotes
Scripts and files for the associated FYP: Comprehensive analysis of stop codon usage in bacteria and archea: "like evidence has suggested in eukaryotes, do prokaryotes select for fail safe stop codons?"

MF_FYP_ASC_Prokaryotes
Scripts and files for the associated FYP "As evidence has suggested in eukaryotes, do prokaryotes select for additional fail-safe stop codons?"
This repository contains the source code for running the analysis described in the project.
Instructions:
• All files ran from root directory
• Scripts are labelled numerically in order of methodology (1 to 26)
EMBL files:
• archea_filtered_genomes
• bacterial_filtered_genomes
Python/R scripts:
A spreadsheet with more detailed script information (Python and R) can be found in the document 'scripts_explained'.
Python:
1.	Extracting untranslated regions (UTRs) of coding sequence (CDS)
2.	Calculating actual observed ASC frequencies in 3’UTR of CDS
3.	Tag null rates simulation
4.	Tga null rates simulation
5.	Taa null rates simulation
6.	Calculating GC3 content
7.	Extracting UTRs of non-coding (NCDS)
8.	Calculating expected ASC frequencies in 3’UTR of NCDS (null 2)
9.	Extracting UTRs of out-of-frame (OOF) sequence
10.	Calculating expected ASC frequencies in 3’UTR of OOF (null 3)
11.	Calculating dinucleotide content
12.	Producing genome master table
13.	Merging files for R analysis
14.	Calculating actual observed TGT frequencies in 3'UTR of CDS
15.	Calculating actual observed TGG frequencies in 3'UTR of CDS
16.	Calculating expected TGT frequencies in 3'UTR of NCDS
17.	Calculating expected TGG frequencies in 3'UTR of NCDS
18.	Calculating expected TGT frequencies in 3'UTR of OOF seq
19.	Calculating expected TGG frequencies in 3'UTR of OOF seq
20.	Merges csv files (TGT and TGG related) for R analysis
R:
14.	Statistical analysis (two-tailed binomial test) and representation of simulated seq null 1 (z scores) histogram (figure_1)
15.	Statistical analysis (Spearman's rank correlation) and representation GC3 content with z scores scatter graph (figure_2)
16.	Statistical analysis (paired t test) and representation of NCDS null 2 (observed - expected) histogram (figure_3)
17.	Statistical analysis (paired t test) and representation of OOF null 3 (observed - expected) histogram (figure_4)
18.	Statistical analysis (paired t test) and representation of dinuclotide effects (observed vs expected) histogram (figure_5)
19.	Statistical analysis (paired t test) and representation of observed TGT and TGG rates vs null 2 and 3 (NCDS and OOF) histogram (appendix 2)
FASTA folders/files explained:
The working repository contains 3 folders of FASTA files:
• UTR_CDS: This folder contains FASTA files for filtered bacteria dataset. Sequences contained in the FASTA files are the 3’ UTR of protein coding sequences with information for all qualifying genes that satisfy filters.
• UTR_NCDS: This folder contains FASTA files for filtered bacteria dataset. Sequences contained in the FASTA files are the 3’ UTR of non-protein coding sequences with information for all qualifying genes that satisfy filters.
• UTR_OOF: This folder contains FASTA files for filtered bacteria dataset. Sequences contained in the FASTA files are the 3’ UTR of non-protein coding sequences with information for all qualifying genes that satisfy filters.
