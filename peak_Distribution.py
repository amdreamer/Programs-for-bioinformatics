#!/usr/bin/env python
# Author: Meng-Ran Wang <wangmengran@picb.ac.cn>
# Info: Chip-seq, CHIRP-seq or RIP CHART-seq are able to enrichment the regionthis script is used to calculate the peak distribution in genes, whether the peak is near 3'UTR or 5'UTR.
# Usage: ./peak_Distribution.py list.txt output_prefix
# input file format: columns in list.txt ('GeneSymbol','chr.peak','Start.peak','End.peak','Length.peak','FoldChange.peak','-lgP.value','GeneID','chr.gene','strand','Start.trans','End.trans','Start.cds','End.cds','Exon.num','Start.exon','End.exon')
from __future__ import division
import re
import sys,os
import pandas as pd
import numpy as np 

in_file = sys.argv[1]
prefix = sys.argv[2]
print in_file

def _judge_region(Strans,Etrans,Scds,Ecds,Strand,SEUTR5,EEUTR5,SEcds,EEcds,SEUTR3,EEUTR3,LenUTR5,LenCds,LenUTR3,pos):
    """To judge the genomic posion is belong to 5'UTR, CDS or 3'UTR and calculate the relative position of the genomics position """
    RelPos = 0
    if Strand == '+':
        if pos < Scds: #"""To judge the genomic posion is belong to 5'UTR"""
        	tag = '5UTR'
        	for i in range(0,SEUTR5.size):
        		if pos > SEUTR5[i] and pos <= EEUTR5[i]:
        			EEUTR5[i] = pos
        			RelPos = format(sum(EEUTR5[0:i+1]-SEUTR5[0:i+1])/LenUTR5,'.4f')
        			break
        else:
            if pos <= Ecds: #"""To judge the genomic posion is belong to CDS"""
                tag = 'CDS'
                for i in range(0,SEcds.size):
                    if pos > SEcds[i] and pos <= EEcds[i]:
                        EEcds[i] = pos; 
                        RelPos = format(sum(EEcds[0:i+1]-SEcds[0:i+1])/LenCds,'.4f')
                        break
                
            else: #"""To judge the genomic posion is belong to 3'UTR"""
            	tag = '3UTR'
                for i in range(0,SEUTR3.size):
                	if pos > SEUTR3[i] and pos <= EEUTR3[i]:
	                	EEUTR3[i] = pos;
	                	RelPos = format(sum(EEUTR3[0:i+1]-SEUTR3[0:i+1])/LenUTR3,'.4f')
	                	break
	            
    if Strand == '-':
        if pos > Ecds: #"""To judge the genomic posion is belong to 5'UTR"""
        	tag = '5UTR'
        	for i in range(0,SEUTR5.size):
        		if pos > SEUTR5[i] and pos <= EEUTR5[i]:
        			SEUTR5[i] = pos
        			RelPos = format(sum(EEUTR5[i:SEUTR5.size+1]-SEUTR5[i:SEUTR5.size+1])/LenUTR5,'.4f')
        			break
        else:
            if pos >= Scds: #"""To judge the genomic posion is belong to CDS"""
            	tag = 'CDS'
                for i in range(0,SEcds.size):
                    if pos > SEcds[i] and pos <= EEcds[i]:
                        SEcds[i] = pos; 
                        RelPos = format(sum(EEcds[i:SEcds.size+1]-SEcds[i:SEcds.size+1])/LenCds,'.4f')
                        break
            else: #"""To judge the genomic posion is belong to 3'UTR"""
            	tag = '3UTR'
                for i in range(0,SEUTR3.size):
                	if pos > SEUTR3[i] and pos <= EEUTR3[i]:
	                	SEUTR3[i] = pos;
	                	RelPos = format(sum(EEUTR3[i:SEUTR3.size+1]-SEUTR3[i:SEUTR3.size+1])/LenUTR3,'.4f')
	                	break
    return tag, RelPos 

def _calculate_region_len(Strans,Etrans,Scds,Ecds,Sexon,Eexon,Strand):
    """To calculate the 5UTR CDS and 3UTR length for each gene and in average"""
    SEUTR5= EEUTR5= SEcds= EEcds= SEUTR3= EEUTR3 = np.array([0])
    LenUTR5= LenCds= LenUTR3 = 0
    if Strand == '+':
        SE, EE = np.array(Sexon[:-1].split(","),dtype=np.int32), np.array(Eexon[:-1].split(","),dtype=np.int32)
        # replace the start of exon with start of CDS, end of exon with end of CDS, at proper exon position!!
        for i in range(0,SE.size):
            if Scds > SE[i] and Scds <= EE[i]:
            	EEUTR5 = EE[0:i+1]; EEUTR5[i] = Scds; SEUTR5 = SE[0:i+1]; LenUTR5 = sum(EEUTR5-SEUTR5)
                SE[i] = Scds; SEcds = np.delete(SE,range(0,i)); EEcds= np.delete(EE,range(0,i))
                break
        for i in list(reversed(range(0,SE.size))):
            if Ecds <= EE[i] and Ecds > SE[i]:
            	SEUTR3 = SE[i:SE.size+1]; SEUTR3[0] = Ecds; EEUTR3 = EE[i:SE.size+1]; LenUTR3 = sum(EEUTR3-SEUTR3)
                EE[i] = Ecds; EEcds = np.delete(EE,range(i+1,SE.size)); SEcds = np.delete(SE,range(i+1,SE.size))
                break
        LenCds = sum(EEcds-SEcds)
    if Strand == '-':
        SE, EE = np.array(Sexon[:-1].split(","),dtype=np.int32), np.array(Eexon[:-1].split(","),dtype=np.int32)
        for i in range(0,SE.size):
            if Scds > SE[i] and Scds <= EE[i]:
            	EEUTR3 = EE[0:i+1]; EEUTR3[i] = Scds; SEUTR3 = SE[0:i+1]; LenUTR3 = sum(EEUTR3-SEUTR3)
                SE[i] = Scds; SEcds = np.delete(SE,range(0,i)); EEcds = np.delete(EE,range(0,i))
                break
        for i in list(reversed(range(0,SE.size))):
            if Ecds <= EE[i] and Ecds > SE[i]:
            	SEUTR5 = SE[i:SE.size+1]; SEUTR5[0] = Ecds; EEUTR5 = EE[i:SE.size+1]; LenUTR5 = sum(EEUTR5-SEUTR5)
                EE[i] = Ecds; EEcds = np.delete(EE,range(i+1,SE.size)); SEcds = np.delete(SE,range(i+1,SE.size))
                break
        LenCds = sum(EEcds-SEcds)
    return SEUTR5, EEUTR5, SEcds, EEcds, SEUTR3, EEUTR3, LenUTR5, LenCds, LenUTR3

# read in data.frame
data = pd.read_table(in_file, header = 0, sep = '\s+', names=['GeneSymbol','chr.peak','S.peak','E.peak','Len.peak','FC.peak','-lgP.value','GeneID','chr.gene','strand','S.trans','E.trans','S.cds','E.cds','Exon.num','S.exon','E.exon'])
# establish the output table
UTR5tab, UTR3tab, CDStab = pd.Series(), pd.Series(), pd.Series()
Len_tab = pd.DataFrame()
# traverse every line
for index, row in data.iterrows():
    Strans,Etrans,Scds,Ecds,Sexon,Eexon,Strand,GeneSymbol = row['S.trans'],row['E.trans'],row['S.cds'],row['E.cds'],row['S.exon'],row['E.exon'],row['strand'],row['GeneSymbol']
    if row['chr.peak'] == row['chr.gene']:
	    SEUTR5, EEUTR5, SEcds, EEcds, SEUTR3, EEUTR3, LenUTR5, LenCds, LenUTR3 = _calculate_region_len(Strans,Etrans,Scds,Ecds,Sexon,Eexon,Strand)
	    Len_tab = Len_tab.append(pd.DataFrame(data={'LenUTR5':LenUTR5,'LenCds':LenCds,'LenUTR3':LenUTR3},index=[str(GeneSymbol)]))
	    # transverse every position
	    for pos in range(row['S.peak'],row['E.peak']+1):
	        tag, RelPos = _judge_region(Strans,Etrans,Scds,Ecds,Strand,SEUTR5,EEUTR5,SEcds,EEcds,SEUTR3,EEUTR3,LenUTR5,LenCds,LenUTR3,pos)
	        if tag == '5UTR':
	            UTR5tab = UTR5tab.append(pd.Series([RelPos],index=[str(GeneSymbol)],dtype='float64'))
	            continue
	        if tag == 'CDS':
	            CDStab = CDStab.append(pd.Series([RelPos],index=[str(GeneSymbol)],dtype='float64'))
	            continue
	        if tag == '3UTR':
	            UTR3tab = UTR3tab.append(pd.Series([RelPos],index=[str(GeneSymbol)],dtype='float64'))
Len_tab = Len_tab[~Len_tab.index.duplicated()]
meanLen = Len_tab.mean(axis = 0,numeric_only = True)

meanLen.to_csv(str(prefix)+'mean_UTR5UTR3CDS_len.csv')
Len_tab.to_csv(str(prefix)+'gene_UTR5UTR3CDS_len.csv')
UTR5tab.to_csv(str(prefix)+'UTR5_relative_pos.csv')
UTR3tab.to_csv(str(prefix)+'UTR3_relative_pos.csv')
CDStab.to_csv(str(prefix)+'CDS_relative_pos.csv')  


# a total table used for plotting the peak distribution of all UTR5 UTR3 and cds region. 
print meanLen
UTR5tab = UTR5tab * meanLen['LenUTR5']
CDStab = meanLen['LenUTR5'] + CDStab * meanLen['LenCds']
UTR3tab = meanLen['LenUTR5'] + meanLen['LenCds'] + UTR3tab * meanLen['LenUTR3'] 
total = UTR5tab.append(CDStab).append(UTR3tab)
total.to_csv(str(prefix)+'peakDistribution.csv')
