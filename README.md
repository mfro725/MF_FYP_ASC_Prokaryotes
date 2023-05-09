# MF_FYP_ASC_Prokaryotes

Scripts and files for the associated FYP carrying out comprehensive analysis of stop codon usage in 728 bacterial and archaic genomes: "Comprehensive analysis of additional stop codon usage in bacteria and archaea: an error-proofing role is minimal, but it can depend on the analysis"


This repository contains the source code for running the analysis described in the project.

====================================================================================

**Instructions:**

• All files ran from root directory

• Scripts are labelled numerically in chronological order of methodology conducted (1 to 26)


**EMBL files:**

• archea_filtered_genomes

• bacterial_filtered_genomes

**Python/R scripts:**

A spreadsheet with more detailed script information (Python and R) can be found in the document 'scripts_explained'.

**Python:**

1.	Extracting untranslated regions (UTRs) of coding sequence (CDS)
3.	Calculating actual observed ASC frequencies in 3’UTR of CDS
4.	Tag null rates simulation
5.	Tga null rates simulation
6.	Taa null rates simulation
7.	Calculating GC3 content
8.	Extracting UTRs of non-coding (NCDS)
9.	Calculating expected ASC frequencies in 3’UTR of NCDS (null 2)
10.	Extracting UTRs of out-of-frame (OOF) sequence
11.	Calculating expected ASC frequencies in 3’UTR of OOF (null 3)
12.	Calculating dinucleotide content
13.	Producing genome master table
14.	Merging files for R analysis
15. Calculating actual observed TGT frequencies in 3'UTR of CDS
16. Calculating actual observed TGG frequencies in 3'UTR of CDS
17.	Calculating expected TGT frequencies in 3'UTR of NCDS
18.	Calculating expected TGG frequencies in 3'UTR of NCDS
19.	Calculating expected TGT frequencies in 3'UTR of OOF seq
20.	Calculating expected TGG frequencies in 3'UTR of OOF seq
21.	Merges csv files (TGT and TGG related) for R analysis


**R:**

21.	Statistical analysis (two-tailed binomial test) and representation of simulated seq null 1 (z scores) histogram (figure_1)
22.	Statistical analysis (Spearman's rank correlation) and representation GC3 content with z scores scatter graph (figure_2)
23.	Statistical analysis (paired t test) and representation of NCDS null 2 (observed - expected) histogram (figure_3)
24.	Statistical analysis (paired t test) and representation of OOF null 3 (observed - expected) histogram (figure_4)
25.	Statistical analysis (paired t test) and representation of dinuclotide effects (observed vs expected) histogram (figure_5)
26.	Statistical analysis (paired t test) and representation of observed TGT and TGG rates vs null 2 and 3 (NCDS and OOF) histogram (appendix 2)

**FASTA folders/files explained:**

The code in this working repository generates 3 folders of FASTA files:

• UTR_CDS: This folder contains FASTA files for filtered bacteria dataset. Sequences contained in the FASTA files are the 3’ UTR of protein coding sequences with information for all qualifying genes that satisfy filters.

• UTR_NCDS: This folder contains FASTA files for filtered bacteria dataset. Sequences contained in the FASTA files are the 3’ UTR of non-protein coding sequences with information for all qualifying genes that satisfy filters.

• UTR_OOF: This folder contains FASTA files for filtered bacteria dataset. Sequences contained in the FASTA files are the 3’ UTR of non-protein coding sequences with information for all qualifying genes that satisfy filters.

