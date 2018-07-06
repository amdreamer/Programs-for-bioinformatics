# Programs-for-bioinformatics
-----------------------------
grep_Lines.py

usage: grepGenes.py [-h] [-f F] searchFile matchFile outputFile

Greping lines from a file (search file) based on certain gene names that provided in another file (match file). 
Please set the Nth field you want to match by -f, ex: -f 2.

positional arguments:
  searchFile
  matchFile
  outputFile

optional arguments:
  -h, --help  show this help message and exit
  -f F        set which field in match file used for matching the search file.
  
 example:
 Grep lines form refFlat files(gene name is in the first field) based on genes in genelist.txt (one gene per line)
 grepGenes.py genelist.txt /picb/rnomics1/database/Human/hg38/refFlat.txt output.txt -f 0
