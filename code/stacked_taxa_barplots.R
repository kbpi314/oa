library(qiime2R)
library(phyloseq)
library(dplyr)
library(ggplot2)

### Themes and colors ###

# background theme
bkg <- theme_q2r() +
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 14, face = "bold", color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0,4,0,0), "mm"))) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 7)) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.spacing.x=unit(1, "lines")) +
  theme(strip.text = element_text(size = 8, face = "bold", color = "black")) +
  theme(plot.title = element_text(size=14, color= "black", face ="bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

col1 <- c("#DDDDDD", "#C4C4C4", "#acacac","#909090", "#757575","#191919","#720510","#a40818","#bf2131", "#C62702","#DD2C03","#FE7926", "#fdb529", "#F3CE5A","#F8DE7E","#394833", "#196919","#528652","#8caf7c", "#b6d7a8", "#bc7e9e","#A55780","#8f2e61","#870047","#6E1945", "#022f59","#005496", "#1075b7","#6ea1ca", "#b3cde0", "#D1E1EC")
col4 <- c("#DDDDDD", "#bc7e9e", "#915d78","#870047","#5E0031", "#C62702", "#022f59","#005496","#6ea1ca", "#b3cde0", "#528652")
col3 <- c("#DDDDDD", "#191919", "#bc7e9e", "#915d78","#870047","#5E0031", "#C62702", "#FE7926", "#F3CE5A", "#022f59","#005496","#6ea1ca", "#b3cde0", "#193843", "#528652")


metadata <- read_q2metadata("qiime_mapping_file_2024Feb8.tsv")
table <- read_qza("taxa_table_L7.qza")$data
taxonomy <- read_qza("taxonomy.qza")$data %>% parse_taxonomy()  

taxasums <- summarize_taxa(table, taxonomy)

p <- taxa_barplot(table, metadata, ntoplot=30) +
  ylab("Relative abundance (%)") + scale_fill_manual(values = col1) + bkg
print(p)
ggsave("metaphlan_taxa_plot_phylum_sorted_ParticipantType.png", height=2200, width=5000, units = 'px') 


taxa_barplot<-function(features, metadata, category, normalize, ntoplot, sort_by="none"){
  library(tidyverse)
  
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
  
  if(!missing(category)){bplot<-bplot + facet_grid(~get(category), scales="free_x", space="free", labeller = as_labeller(wrap_text))}
  
  return(bplot)
}
