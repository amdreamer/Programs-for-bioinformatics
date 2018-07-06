#!/bin/bash
# 2018.05.10 17:40
# One gene may has more than one isoforms. This program was wrote for greping longest isoform for each gene in a fasta file.

Usage() {
echo"Usage: $0 $2"
}
# judge file's format

# read fasta file
file=$1
basename=$(basename $file .fa)
seqkit fx2tab $file > $basename.tab
#echo "$basename.tab"
#head $basename.tab
# header info: which field is geneName: default $1
awk -F'[|"\t"]' '{print $1"\t"$2}' $basename.tab|sort > $basename.2.tab
#echo "$basename.2.tab"
#grep "AAGAB" $basename.2.tab
# find the longest isofrom based on the sequence length.
awk '{if (len[$1] < length($2)){len[$1]=length($2);a[$1]=$2;}}END{for(i in a){print i"\t"a[i];}}' $basename.2.tab > $basename.longest.tab
#echo "$basename.longest.tab"
#grep "AAGAB" $basename.longest.tab
# transfer into fasta format (.fa)
seqkit tab2fx $basename.longest.tab > $basename.longest.fa
#echo "$basename.longest.fa"
#head $basename.longest.fa
# remove tmp files
rm $basename.tab $basename.2.tab
