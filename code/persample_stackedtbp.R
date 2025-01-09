# stacked taxabarplot from AC
library(tidyverse)
library(ggplot2)

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
bkg <- theme_q2r() +
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) + 
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 14, face = "bold", color = "black")) +
  # theme(axis.title.y = element_blank()) + #element_text(margin = unit(c(0,4,0,0), "mm"))) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 7)) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.spacing.x=unit(1, "lines")) +
  theme(strip.text = element_text(size = 8, face = "bold", color = "black")) +
  theme(strip.text = element_blank()) 
  # theme(plot.title = element_text(size=14, color= "black", face ="bold")) +
  # theme(plot.title = element_text(hjust = 0.5)) # +
#  theme(plot.title = element_blank()) 
# theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

col3 <- c("#DDDDDD", "#191919", "#bc7e9e", "#915d78","#870047","#5E0031", "#C62702", "#FE7926", "#F3CE5A", "#022f59","#005496","#6ea1ca", "#b3cde0", "#193843", "#528652")
col4 <- c("#DDDDDD", "#bc7e9e", "#915d78","#870047","#5E0031", "#C62702", "#022f59","#005496","#6ea1ca", "#b3cde0", "#528652")
col1 <- c("#DDDDDD", "#C4C4C4", "#acacac","#909090", "#757575","#191919","#720510","#a40818","#bf2131", "#C62702","#DD2C03","#FE7926", "#fdb529", "#F3CE5A","#F8DE7E","#394833", "#196919","#528652","#8caf7c", "#b6d7a8", "#bc7e9e","#A55780","#8f2e61","#870047","#6E1945", "#022f59","#005496", "#1075b7","#6ea1ca", "#b3cde0", "#D1E1EC")
col2 <- c("#DDDDDD", "#acacac", "#757575","#191919","#720510","#a40818","#DD2C03","#FE7926", "#fdb529", "#F3CE5A","#F8DE7E","#394833", "#196919","#8caf7c", "#b6d7a8", "#bc7e9e","#870047","#6E1945", "#022f59","#005496", "#1075b7","#6ea1ca", "#b3cde0")

taxa_barplot<-function(features, metadata, category, normalize, ntoplot, col, sort_by="none"){
  q2r_palette<-c(
    "blue4",
    "olivedrab",
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "aquamarine",
    "yellow3",
    "coral",
    "grey"
  )
  q2r_palette <- col
  
  if(missing(ntoplot) & nrow(features)>10){ntoplot=10} else if (missing(ntoplot)){ntoplot=nrow(features)}
  if(missing(normalize)){normalize<-"percent"}
  if(normalize=="percent"){features<-make_percent(features)} else if(normalize=="proportion"){features<-make_proportion(features)}
  if(missing(metadata)){metadata<-data.frame(SampleID=colnames(features))}
  if(!"SampleID" %in% colnames(metadata)){metadata <- metadata %>% rownames_to_column("SampleID")}
  if(!missing(category)){
    if(!category %in% colnames(metadata)){message(stop(category, " not found as column in metdata"))}
  }
  
  plotfeats<-names(sort(rowMeans(features), decreasing = TRUE)[1:ntoplot]) # extract the top N most abundant features on average
  
  suppressMessages(
    suppressWarnings(
      fplot<-
        features %>%
        as.data.frame() %>%
        rownames_to_column(var="Taxon") %>%
        gather(-Taxon, key="SampleID", value="Abundance") %>%
        mutate(Taxon=if_else(Taxon %in% plotfeats, Taxon, "Remainder")) %>%
        group_by(Taxon, SampleID) %>%
        summarize(Abundance=sum(Abundance)) %>%
        ungroup() %>%
        mutate(Taxon=factor(Taxon, levels=rev(c(plotfeats, "Remainder")))) %>%
        left_join(metadata)
    ))
  bplot<-
    ggplot(fplot, aes(x=reorder(SampleID, ave(setNames(Abundance, Taxon), SampleID, FUN = function(x) replace(x[sort_by], is.na(x[sort_by]), 0))), y=Abundance, fill=Taxon)) +
    geom_bar(stat="identity") +
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_cartesian(expand=FALSE) +
    xlab("Sample") +
    ylab("Abundance")
  if(ntoplot<=10){bplot<-bplot+scale_fill_manual(values=rev(q2r_palette), name="")}
  if(!missing(category)){bplot<-bplot + facet_grid(~get(category), scales="free_x", space="free")}#, labeller = as_labeller(wrap_text))}
  return(bplot)
}

# read in metadata
# metadata=read_delim(file='/Users/KevinBu/Desktop/clemente_lab/Projects/oa/inputs/df_mapdf_meta.tsv', delim='\t')#,row
# metadata=metadata[-1,] # only if using Q2 mapping file 


# plots
#visits = c('saliva_post','saliva_pre','stool_post','stool_pre','18mo')
#levels = c('L6')
#all_pairs <- c(outer(visits, levels, FUN = paste, sep = "_")) 
all_pairs <- c('saliva_post','saliva_pre','stool_post','stool_pre') 
# "1_2wk_L2" "6wk_L2"   "6mo_L2"   "12mo_L2"  "18mo_L2"  "1_2wk_L6" "6wk_L6"   "6mo_L6"   "12mo_L6"  "18mo_L6" 

sortby = c(
  'Streptococcus','Prevotella',
  'Bacteroides_H','Bacteroides_H'
)


for (i in seq(1:length(all_pairs))){
  s=all_pairs[i]
  sp_taxasums = read.csv(file=paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs35/df_q2R_',s,'.tsv'), 
                         sep='\t',
                         row.names='Taxon')
  print(length(sp_taxasums))
  
  nplot=15
  p <- taxa_barplot(sp_taxasums, 
                    #metadata, 
                    #category="Group", 
                    ntoplot=nplot, 
                    col=col2, 
                    sort_by=sortby[i],
                    normalize='none') + 
    ylab('Relative abundance (%)') + 
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = col2[1:(nplot+1)]) +
    bkg +
    ggtitle(paste0('Taxabarplot_',s)) 
    
  
  pdf(file = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs35/stacked_barplot_',s,'.pdf'), height = 4, width = 6)
  plot(p)
  dev.off()
}
  
  

