#!/bin/bash
# 2018.05.10 17:40
# One gene may has more than one isoforms. This program was wrote for greping longest isoform for each gene in a fasta file.
# getLongestIsoform.sh input_format input_file output_file

Usage() {
echo"$0: getLongestIsoform.sh input_format input_file"
ecgo"Input file format support fasta or genePred format so far. correspondingly, filename= *.fa or *.txt"
}
# judge file's format
if $1 == "fasta"; then
  # read fasta file
  file = $2
  basename = $(basename $file .fa)
  seqkit fx2tab $file > $basename.tab
  # header info: which field is geneName: default $1
  awk -F'[|"\t"]' '{print $1"\t"$2}' $basename.tab|sort > $basename.2.tab
  # find the longest isofrom based on the sequence length.
  awk '{if (len[$1] < length($2)){len[$1]=length($2);a[$1]=$2;}}END{for(i in a){print i"\t"a[i];}}' $basename.2.tab > $basename.longest.tab
  # transfer into fasta format (.fa)
  seqkit tab2fx $basename.longest.tab > $basename.longest.fa
  # remove tmp files
  rm $basename.tab $basename.2.tab
fi
if $1 == "genePred"; then
 file = $2
 basename = $(basename $file .txt)
 awk '{if (lenTrans[$1]<=$6-$5 && lenCDS[$1]<=$8-$7) {len[$1]=$6-$5;lenCDS[$1]=$8-$7;a[$1]=$0;}}END{for(i in a){print a[i];}}' $2 > $basename.longest.txt
fi
