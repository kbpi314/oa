################################# 
## R script                    ##
## Project: OA                 ##
## Beta diversity pcoa plots   ##
## Author: KB                  ##
#################################

############################################################################
############################################################################

### Load libraries ###

#library(phyloseq)
#library(vegan)
#library(ade4)
#library(PMCMR)
#library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)

############################################################################
############################################################################
############################################################################

### Beta diversity pcoa ###

# background theme
bkg <- theme_bw() +
  theme(axis.text = element_text(size = 24, color = "black")) +
  theme(axis.title = element_text(size = 32, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black"))+ #, face = "bold")) +
  # theme(legend.title = element_blank()) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black"))
  theme(legend.justification = "right")# +
  #theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))

# function to specify that axis labels have 2 decimal places
f.dec <- function(x){
  format(round(x, 2), nsmall = 2)
}

# directory for storing files
outdir = "/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs32/"

# list of distance methods
dists <- c('bray_curtis')


paired_dirs = c(
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_response/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_noresponse/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_response/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_noresponse/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_response/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_noresponse/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_response/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_noresponse/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_response/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_noresponse/"
)

paired_p = c(
  '0.985', # stool otu
  '0.960', # res
  '0.913', # nr
  
  '0.490', # saliva otu
  '0.680',
  '0.970',
  
  '0.335', # meta stool
  '0.405',
  '0.375',
  
  '0.815', # meta saliva
  '0.900',
  '0.865',
  
  '0.005', # meta plasma
  '0.035',
  '0.115'
  
)

unpaired_p = c(
  '0.580', # saliva otu pre
  '0.330', # post

  '0.950', # stool otu
  '0.150',
  
  '0.415', # saliva meta
  '0.705',
  
  '0.955', # plasma meta
  '0.805',
  
  '0.635', # stool meta
  '0.280'
)

unpaired_dirs = c(
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_pre/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_post/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_pre/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_post/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_pre/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_post/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_pre/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_post/",
  
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_pre/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_post/"
  
)

col1 <- c("#929aab", "#ce2525")
col1unpaired <- c("green4", "lightblue")
col2 <- c("#f3a333", "#0074e4", "#8f8787")

for (i in 1:length(paired_dirs)){
  dir=paired_dirs[i]
  # load data
  df <- read.delim(file=paste0(dir,"pcoa.tsv"),
                   row.names=1)

  # drop unpaired HostSubjectId
  df <- df |>
    filter(n() > 1, .by = HostSubjectId)
  
  # capitalize Timepoint
  df$Timepoint = toTitleCase(df$Timepoint)
  
  # order factors for legend
  df$Timepoint <- factor(df$Timepoint, levels=c('Pre', 'Post'))

  for (j in seq_along(dists)) {
    # create filenames
    filename_plot = paste(basename(dir),"bdiv", dists[j], "plot.pdf", sep = "_")
    
    # plot beta diversity
    pv = paired_p[i]
    p <- ggplot() + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
      geom_point(data = df, aes(x = PC1, y = PC2, color = Timepoint), size=4) + #shape = sib_02),size=4) +
      scale_color_manual(values = c("Pre" = col1[1], "Post" = col1[2])) + #col1, labels = c("Unaffected", "RA")) +
      bkg + guides(shape = "none") + 
      ggtitle(paste0('p=',pv)) +
      scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
      scale_y_continuous(labels = f.dec)   # 2 decimal places on y-axis
    
    # save plot
    fp = paste(outdir, filename_plot, sep = "")
    pdf(file = fp, height = 6, width = 8)
    plot(p)
    dev.off()
  }
  
}


bkg <- theme_bw() +
  theme(axis.text = element_text(size = 24, color = "black")) +
  theme(axis.title = element_text(size = 32, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black"))+ #, face = "bold")) +
  theme(legend.title = element_blank()) +
  # theme(legend.title = element_text(size = 24, face = "bold", color = "black"))
  theme(legend.justification = "right")# +


for (i in 1:length(unpaired_dirs)){
  dir=unpaired_dirs[i]
  # load data
  df <- read.delim(file=paste0(dir,"pcoa.tsv"),
                   row.names=1)
  
  # dropna
  df = df[!is.na(df$WOMAC_P_Response),]
  
  # order factors for legend
  df$WOMAC_P_Response <- factor(df$WOMAC_P_Response, levels = c('Response', 'No response'))
  
  for (j in seq_along(dists)) {
    # create filenames
    filename_plot = paste(basename(dir),"bdiv", dists[j], "plot.pdf", sep = "_")
    
    # plot beta diversity
    pv = unpaired_p[i]
    p <- ggplot() + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
      geom_point(data = df, aes(x = PC1, y = PC2, color = WOMAC_P_Response), size=4) + #shape = sib_02),size=4) +
      scale_color_manual(values = c("Response" = col1unpaired[1], "No response" = col1unpaired[2])) + #col1, labels = c("Unaffected", "RA")) +
      bkg + guides(shape = "none") + 
      ggtitle(paste0('p=',pv)) +
      theme(title = element_text(size = 18, color = "black")) +
      scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
      scale_y_continuous(labels = f.dec)   # 2 decimal places on y-axis
    
    # save plot
    fp = paste(outdir, filename_plot, sep = "")
    pdf(file = fp, height = 6, width = 8)
    plot(p)
    dev.off()
  }
}

