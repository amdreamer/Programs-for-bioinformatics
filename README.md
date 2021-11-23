# Programs-for-bioinformatics
## [grep_Lines.py](./grep_Lines.py) 
### Usage: grepGenes.py [-h] [-f F] searchFile matchFile outputFile

Greping lines from a file (search file) based on certain gene names that provided in another file (match file). 

Please set the Nth field you want to match by -f, ex: -f 2.

**positional arguments:**

  searchFile，matchFile，outputFile

**optional arguments:**

  -h, --help  show this help message and exit
  
  -f F        set which field in match file used for matching the search file.
  
**example:**
 
（1）Grep lines form refFlat files(gene name is in the first field) based on genes in genelist.txt (one gene per line)
 
 grepGenes.py genelist.txt /picb/rnomics1/database/Human/hg38/refFlat.txt output.txt -f 0
 
 ## [grep_longest_isoform.sh](./grep_longest_isoform.sh)
 ### Usage: grep_longest_isoform.sh input_format input_file
 
 One gene may has more than one isoforms. This program was wrote for greping longest isoform for each gene in a fasta file.
 The input format are supporting "fasta" and "genePred" so far.

**example:**

 getLongestIsoform.sh fasta 1564_genelist_cDNA.fa
 getLongestIsoform.sh genePred /Human/hg38/170823_annotation/refFlat.txt

 ## [peak_Distribution.py](./peak_Distribution.py)
 ### Usage: ./peak_dist_test.py list.txt output_prefix

 input file (list.txt) format: 
 
 columns in list.txt ('GeneSymbol', 'chr.peak', 'Start.peak', 'End.peak', 'Length.peak', 'FoldChange.peak', 
 '-lgP.value', 'GeneID', 'chr.gene', 'strand', 'Start.trans', 'End.trans', 'Start.cds', 'End.cds', 'Exon.num', 'Start.exon', 'End.exon' )
 
  ## [sra_download.py](./sra_download.py)
  ### Usage: ./sra_download.py 
  #### -sraID_file SRAID_FILE
  ####                       passing a file with sra IDs, one ID for each line.
  #### -ID ID [ID ...]       Instead of passing a file with sraID for each line,
  ####                       user can also input sraID directly by using this
  ####                       parameter. ex: -ID SRR2121685 SRR2121686 SRR2121687
  This program can be used for batch downloading fastq files from GEO, using sra IDs.
  
  ## [survival_analysis.R] (./survival_analysis.R)
  ### Survival analysis of TCGA patients integrating gene expression (MS,proteome) data. This script is similar with https://www.biostars.org/p/153013/.
  #### to do survival analysis we need three main things:
  #### time: this is the time till an event happens
  #### status: this indicates which patients have to be kept for the analysis
  #### event: this tells i.e. which patients have the gene up- or down-regulated or have no changes in expression
  #### Since we want to do censored analysis, we need to have something to censor the data with. For example, if a patient has no death data BUT there is a date to last followup it means that after that day we know nothing about the patient, therefore after that day it cannot be used for calculations/Kaplan Meier plot anymore, therefore we censor it. so now we need to create vectors for both 'time to new tumor' and 'time to death' that contain also the data from censored individuals.
  
  
  
