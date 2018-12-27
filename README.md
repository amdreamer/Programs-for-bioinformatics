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
