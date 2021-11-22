## calculate survival rate in CPTCA dataset
clinic = read.csv('TCGA_COAD_READ_clinical.tsv',sep="\t")
# replace '-- to NA
clinic[clinic == "\'--"] = NA
colnames(clinic)
clinic_selct=clinic[,c(2,12,26,28)] # select sampleCode, Metastasis or not, alive or death, survival time
colnames(clinic_selct)=c('SampleCode','Metastasis','status','survivalTime')

# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
ind_keep <- grep('days_to_death',colnames(clinic)) # get the column index of days of death
death <- as.matrix(clinic[,ind_keep]) 
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_follow_up',colnames(clinic))
fl <- as.matrix(clinic[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]) {
  if ( sum(is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(death_collapsed,fl_collapsed)
colnames(all_clin) <- c('death_days', 'followUp_days')
# create vector time to death containing values to censor for death
# For example, if a patient has no death data BUT there is a date to last followup it means that after that day we know nothing about the patient, 
# therefore after that day it cannot be used for calculations/Kaplan Meier plot anymore, therefore we censor it. 
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

# create vector for death censoring
table(clinic$vital_status)
# alive dead
# 372   161
all_clin$death_event <- ifelse(clinic$vital_status == 'Alive', 0,1) # alive:0; death:1
#finally add row.names to clinical
all_clin$SampleCode<- clinic$case_submitter_id

# For each genes, select which samples to perform survival test. combining clinic info with protein expression.
## modify the colnames of p_exp to match the sample code.
sampleCode=substring( gsub('\\.','-',colnames(p_exp)),1,12)
colnames(p_exp)=sampleCode
## start survival analysis
percent=vector()
geneName=vector()
for (gene in genelist) {
  #gene='GIPC2'
  #ID='P20290'
  express_gene=data.frame()
  express_gene=p_exp[which(p_exp$GeneSymbol==gene),] 
  express_gene_2=express_gene[,-1]
  express_gene_2=express_gene_2[,order(express_gene_2[1,])]
  express_gene_2[ , colSums(is.na(express_gene_2)) == 0] # rm NAs
  if (nrow(express_gene_2)==0) {
    next;
  }
  # set percentile cut-off for protein expression
  for (percentile in c(0.1,0.25,0.4,0.5)) {
    # percentile=0.1, select top 10%, bottom 10%
    low_SampleCode_10=colnames(express_gene_2[,c(1:round(ncol(express_gene_2)*percentile))])
    low_SampleCode_10=data.frame(cbind(low_SampleCode_10,rep(-1,length(low_SampleCode_10))))
    high_SampleCode_10=colnames(express_gene_2[,c(round(ncol(express_gene_2)*(1-percentile)):ncol(express_gene_2))])
    high_SampleCode_10=data.frame(cbind(high_SampleCode_10,rep(1,length(high_SampleCode_10))))
    colnames(low_SampleCode_10)=colnames(high_SampleCode_10)=c('SampleCode','exp_group') # -1 for low express, 1 for high express
    exp_group=data.frame(rbind(low_SampleCode_10,high_SampleCode_10))
    
    # combine clinic with samples of expression changes
    merge=merge(exp_group,all_clin,by='SampleCode')
    # 2*2 table and chisq.test
    #a=merge[which(merge$exp_group==1 & startsWith(merge$Metastasis,'M')),]
    #b=merge[which(merge$exp_group==1 & startsWith(merge$Metastasis,'noM')),]
    #c=merge[which(merge$exp_group==(-1) & startsWith(merge$Metastasis,'M')),]
    #d=merge[which(merge$exp_group==(-1) & startsWith(merge$Metastasis,'noM')),]
    #table_chisquare=data.frame(c(nrow(a),nrow(c)),c(nrow(b),nrow(d)))
    #colnames(table_chisquare)=c('M','noM')
    #rownames(table_chisquare)=c('highExp','lowExp')
    #write.table(table_chisquare,paste(gene,percentile,'2x2table.txt'),quote = FALSE)
    #p_chisquare=c(p_chisquare,chisq.test(table_chisquare)$p.value)
    #percent=c(percent,percentile)
    #geneName=c(geneName,gene)
    # survival analysis
    fit<- survfit(Surv(new_death, death_event) ~ exp_group, data = merge)
    summary(fit)
    p=ggsurvplot(fit, data = merge,
                 surv.median.line = "hv", # Add medians survival
                 # Change legends: title & labels
                 legend.title = "Expression",legend.labs = c("Low", "High"),
                 # Add p-value and tervals
                 pval = TRUE,pval.size = 5,pval.method=TRUE,
                 # Change censor
                 censor.shape = 124,censor.size = 2,
                 conf.int = FALSE,# 有无置信区间
                 break.x.by = 400, #横轴坐便
                 # Add risk table
                 #risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
                 #palette = c("#E7B800", "#2E9FDF"),
                 palette = c("blue", "purple"),
                 ggtheme = theme_bw(), # Change ggplot2 theme
                 # Change font size, style and color
                 title = paste(gene,"percent",percentile),
                 font.title = c(14, "bold", "black"),
                 font.x = c(14, "bold.italic", "black"),
                 font.y = c(14, "bold.italic", "black"),
                 font.tickslab = c(12, "plain", "black")
    )
    pdf(paste(gene,percentile,"survival.pdf"),width=8,height = 5,onefile=FALSE)
    print(p)
    dev.off()
  }}
