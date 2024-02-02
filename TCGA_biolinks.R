if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("maftools")
library("TCGAbiolinks")
library(tidyverse)
library(maftools)


paad.clin <- GDCquery_clinic("TCGA-COAD")
clin <- paad.clin %>% 
  dplyr::select(c(
    "submitter_id", "days_to_last_follow_up", "vital_status",
    "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
    "ajcc_pathologic_stage", "gender"
  )) %>%
  `names<-`(c("Tumor_Sample_Barcode", "time", "status", "T", "N", "M", "stage", "gender"))
GDCquery_Maf(
  tumor = "COAD",
  save.csv = TRUE,
  directory = "~/Downloads/TCGA",
  pipelines = "mutect2"
)

paad.maf <- read.maf(
  "~/Downloads/TCGA/TCGA-COAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/03652df4-6090-4f5a-a2ff-ee28a37f9301/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz",
  clinicalData = clin,
  isTCGA = TRUE
)

vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = paad.maf, colors = vc_cols, top = 15)
