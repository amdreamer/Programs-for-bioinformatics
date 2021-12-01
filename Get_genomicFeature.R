# check the features of genes, length 3UTR, 5UTR
library('GenomicFeatures')
BiocManager::install("org.Hs.eg.db")
install.packages("~/Downloads/gencode.v37.basic.annotation_ProteinCoding_201.gtf", repos = NULL, type = "source")
library(dplyr)
library(stringr)

# extract three UTR from GENCODE proteincoding genes (201 main transcripts)
TxDb <- makeTxDbFromGFF("~/Downloads/gencode.v37.basic.annotation_ProteinCoding_201.gtf")
#TxDb <- makeTxDbFromGFF("~/Downloads/gencode.v37.basic.annotation_DDR1_223.gtf")
transcripts.data <- transcripts(TxDb, columns=c("tx_name", "gene_id"))
anyDuplicated(elementMetadata(transcripts.data)$tx_name)


tx2gene <- unlist(elementMetadata(transcripts.data)$gene_id)
names(tx2gene) <- elementMetadata(transcripts.data)$tx_name
## for 3UTR
threeUTRs.data <- threeUTRsByTranscript(TxDb, use.names=T)
# mcols(threeUTRs.data)$gene_id <- drop(transcripts.data$gene_id[match(names(threeUTRs.data),  transcripts.data$gene_id)])
names(threeUTRs.data) <- tx2gene[names(threeUTRs.data)]
length_threeUTRs   <- width(ranges(threeUTRs.data))
the_lengths        <- as.data.frame(length_threeUTRs)
# the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
the_lengths=the_lengths[,-1]
the_lengths$geneID=as.character(lapply(strsplit(the_lengths$group_name,'[.]'),`[[`,1))
the_lengths=the_lengths[,-1]
colnames(the_lengths) <- c("threeUTR_Length","ensembl_gene_id")

## for 5UTR
fiveUTRs.data <- fiveUTRsByTranscript(TxDb, use.names=T)
# mcols(threeUTRs.data)$gene_id <- drop(transcripts.data$gene_id[match(names(threeUTRs.data),  transcripts.data$gene_id)])
names(fiveUTRs.data) <- tx2gene[names(fiveUTRs.data)]
length_fiveUTRs   <- width(ranges(fiveUTRs.data))
the_lengths_fiveUTRs        <- as.data.frame(length_fiveUTRs)
# the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
the_lengths_fiveUTRs=the_lengths_fiveUTRs[,-1]
the_lengths_fiveUTRs$geneID=as.character(lapply(strsplit(the_lengths_fiveUTRs$group_name,'[.]'),`[[`,1))
the_lengths_fiveUTRs=the_lengths_fiveUTRs[,-1]
colnames(the_lengths_fiveUTRs) <- c("fiveUTR_Length","ensembl_gene_id")

## map ensembl ID to gene name
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# if gene ID is formated as ENSG*.*, keep content before ".":
genes <- lapply(strsplit(the_lengths$geneID,'[.]'), `[[`, 1)
genes=as.character(genes)
# if gene ID is formated as ENSG*
genes=as.character(the_lengths$ensembl_gene_id)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

# merge gene symbol with UTR length
merge1= merge(the_lengths,G_list,by='ensembl_gene_id')
merge2= merge(merge1,the_lengths_fiveUTRs,by='ensembl_gene_id')
# select the first record for dulicated rows with same gene names
merge2=merge2[!duplicated(merge2$ensembl_gene_id),]
write.table(merge2,'hg38_gene_UTRlength.txt',row.names=FALSE,quote = FALSE,sep = "\t")

# load MS data and plot UTR lentgh by group
df_KD_MS=read.csv('/Users/amdreamer/Documents/博后期间课题/EMT/DDX6KO_MS/diffenence_DDX6KDRNA_DDX6KOprotein.csv')
df_KD_MS[which(df_KD_MS$log2FC.protein<(-1)),]
df_KD_MS$tag3=as.numeric(df_KD_MS$log2FC.protein>1)*2+as.numeric(df_KD_MS$log2FC.protein>(-1) & df_KD_MS$log2FC.protein<1) # up:2, n=15; unchange:1,n=585; down:0, n=53
merge3=merge(df_KD_MS,merge2,by.x='Gene.names',by.y='hgnc_symbol')
ggplot(merge3,aes(x=factor(tag3),y=log2(threeUTR_Length)))+geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+ylim(0,15)
## stat test
wilcox.test(merge3[which(merge3$tag3==0),'threeUTR_Length'],merge3[which(merge3$tag3==1),'threeUTR_Length']) # p-value = 0.5011
wilcox.test(merge3[which(merge3$tag3==0),'threeUTR_Length'],merge3[which(merge3$tag3==2),'threeUTR_Length']) # p-value = 0.4639
wilcox.test(merge3[which(merge3$tag3==1),'threeUTR_Length'],merge3[which(merge3$tag3==2),'threeUTR_Length']) # p-value = 0.9839

ggplot(merge3,aes(x=factor(tag3),y=log2(fiveUTR_Length)))+geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+ylim(0,11)
## stat test
wilcox.test(merge3[which(merge3$tag3==0),'fiveUTR_Length'],merge3[which(merge3$tag3==1),'fiveUTR_Length']) # p-value = 0.734
wilcox.test(merge3[which(merge3$tag3==0),'fiveUTR_Length'],merge3[which(merge3$tag3==2),'fiveUTR_Length']) # p-value = 0.2701
wilcox.test(merge3[which(merge3$tag3==1),'fiveUTR_Length'],merge3[which(merge3$tag3==2),'fiveUTR_Length']) # p-value = 0.02475


# extract three UTR from GENCODE proteincoding genes (201 main transcripts)
TxDb <- makeTxDbFromGFF("~/Downloads/gencode.v37.basic.annotation_201.gtf")
TxDb <- makeTxDbFromGFF("~/Downloads/gencode.v37.basic.annotation_DDR1_223.gtf")
TxDb <- makeTxDbFromGFF("~/Downloads/gencode.v37.basic.annotation.gtf")
transcripts.data <- transcripts(TxDb, columns=c("tx_name", "gene_id"))
anyDuplicated(elementMetadata(transcripts.data)$tx_name)

tx2gene <- unlist(elementMetadata(transcripts.data)$gene_id)
names(tx2gene) <- elementMetadata(transcripts.data)$tx_name

# create SAF-format annotation file from GTF for featurecounts
create_df_names <- function(gr, filename, record.names=NULL){
  df <- data.frame(GeneID= gsub('\\.[0-9]+', '', record.names),
                   Chr=seqnames(gr),
                   Start=start(gr)-1,
                   End=end(gr),
                   # scores=c(rep('.', length(gr))),
                   Strand=strand(gr))
  df <- arrange(df, Chr, Start, End, Strand)
  #colnames(df) = c('GeneID', 'Chr', 'Start', 'End', 'Strand')
  df <- unique(df)
  write.table(df, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
  return(df)
}

## for 3UTR
threeUTRs.data <- unlist(threeUTRsByTranscript(TxDb, use.names=T))
threeUTRs.df = create_df_names(threeUTRs.data, '~/Downloads/gencode.v37.basic.annotation.threeUTRs.SAF', names(threeUTRs.data))
## for 5UTR
fiveUTRs.data <- unlist(fiveUTRsByTranscript(TxDb, use.names=T))
fiveUTRs.df = create_df_names(fiveUTRs.data, '~/Downloads/gencode.v37.basic.annotation.fiveUTRs.SAF', names(fiveUTRs.data))
# for CDS 
CDS.data = unlist(cdsBy(TxDb, use.names=T))
CDS.df=create_df_names(CDS.data, '~/Downloads/gencode.v37.basic.annotation.CDS.SAF', names(CDS.data))

# create SAF-format annotation file from GTF for featurecounts
create_df_names <- function(gr, filename, record.names=NULL){
  df <- data.frame(GeneID= gsub('\\.[0-9]+', '', record.names),
                   Chr=seqnames(gr),
                   Start=start(gr)-1,
                   End=end(gr),
                   # scores=c(rep('.', length(gr))),
                   Strand=strand(gr))
  df <- arrange(df, Chr, Start, End, Strand)
  #colnames(df) = c('GeneID', 'Chr', 'Start', 'End', 'Strand')
  df <- unique(df)
  write.table(df, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
  return(df)
}
