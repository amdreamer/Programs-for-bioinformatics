# plot boxplot with lines connectted paired points
theme_set(theme_bw(16))
  ggplot(df,aes(factor(tag),log2FC)) +
  geom_boxplot() +
  geom_point(alpha=0.3)+ 
  #geom_text(aes(label=ifelse(log2FC>1.2,as.character(Gene.names),'')),hjust=1, vjust=-0.7,position = position_dodge(0.5))+
  geom_line(aes(group=paired,color=factor(trend)),alpha=0.3) +
  theme(legend.position = "none")+
  scale_color_manual(values = c('gray','darkblue','red'))
