################################# 
## R script                    ##
## Alpha diversity             ##
## Author: KB.                 ##
## Date: 12/29/2024            ##
## Last Updated: 12/29/2024    ##
#################################

### Load libraries ###
library(reshape2)
library(phyloseq)
library(vegan)
library(ade4)
library(PMCMR)
library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(cowplot)
# library(ggpubr) # needs car library

############################################################################
############################################################################
############################################################################

### Statistics functions ###

# All plot statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

# Boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

# Whiskers statistics: median, 10th percentile, 90th percentile
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

# Outliers
min.outlier <- function(x) {
  subset(x, quantile(x, prob = c(0.10)) > x)
}

max.outlier <- function(x) {
  subset(x, quantile(x, prob = c(0.90)) < x)
}

############################################################################
############################################################################
############################################################################

# for each comparison
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

paired_cols = list(
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh/",
  c('red2','pink','lightgray','darkgray'), # good
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_response/",
  c('red2','pink','lightgray','darkgray'), # 
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_noresponse/",
  c('pink','lightgray'),# 
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh/",
  c('pink','lightgray','darkgray'), #
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_response/",
  c('red2','pink','lightgray','darkgray'),#
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_noresponse/",
  c('pink','lightgray'),#
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh/",
  c('pink','lightgray','darkgray'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_response/",
  c('red2','pink','lightgray','darkgray'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_noresponse/",
  c('pink','lightgray','darkgray'),
  
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh/",
  c('red2','pink','lightgray','darkgray'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_response/",
  c('pink','lightgray','darkgray'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_noresponse/",
  c('pink','lightgray'),
  # 
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh/",
  c('red4','red2','pink','lightgray','darkgray','black'), # Meta_plasma_adh
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_response/",
  c('red2','pink','lightgray','darkgray'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_noresponse/"
  c('red2','pink','lightgray')
  # 
  # c("black","darkgray","lightgray",'pink','red2','red4')
)

unpaired_cols = list(
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_pre/",
  c('lightgreen','lightblue'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_post/",
  c('lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_pre/",
  c('green4','lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_post/",
  c('green4','lightgreen','lightblue','blue2'),
  
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_pre/",
  c('green4','lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_post/",
  c('lightgreen','lightblue'),
  #
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_pre/",
  c('lightgreen','lightblue'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_post/",
  c('green4','lightgreen','lightblue','blue2'),
  # 
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_pre/",
  c('green4','lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_post/"
  c('green4','lightgreen','lightblue','blue2')
)
diff_dirs = c(
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_diff/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_diff/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_diff/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_diff/",
  "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_diff/"
)
diff_cols = list(
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_stool_adh_diff/",
  c('green4','lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_stool_adh_diff/",
  c('green4','lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_plasma_adh_diff/",
  c('lightgreen','lightblue'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Qiime2_saliva_adh_diff/",
  c('green4','lightgreen','lightblue','blue2'),
  # "~/Desktop/clemente_lab/Projects/oa/outputs/Meta_saliva_adh_diff/"
  c('lightgreen','lightblue')
)

############################################################################
############################################################################
############################################################################

abkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 18, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 12))+#, color = "black")) +
  theme(axis.title.y = element_text(size = 18, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black")) +
  theme(legend.position='none') + 
  theme(title = element_text(size = 12, color = "black")) 


bbkg <- theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 12, color = "black"))+ #, face = "bold")) +
  theme(legend.title = element_blank()) +
  # theme(legend.title = element_text(size = 24, face = "bold", color = "black"))
  theme(legend.justification = "right") + # +
  # theme(legend.position='none') + 
  theme(title = element_text(size = 12, color = "black")) 
  
test = list()
for (i in 1:length(paired_dirs)){
  # alpha plots
  # background theme

  # choose colors
  col1 <- c("#929aab", "#ce2525")
  # col1unpaired <- c("green4", "purple")
  col1unpaired <- c("green4", "lightblue")
  col2 <- c("#f3a333", "#0074e4", "#8f8787")
  
  # choose line types
  line1 <- c("solid", "dashed", "dotted")
  
  # set output directory
  outdir="~/Desktop/clemente_lab/Projects/oa/outputs/jobs31/"
  
  
  dir = paired_dirs[i]
  print(dir)
  
  df_alpha = read.table(paste0(dir,'metadata.tsv'), 
                        sep = '\t', 
                        header = TRUE, 
                        row.names = 1,
                        check.names = FALSE, 
                        na.strings = "NA")
  
  #drop na
  df_alpha = df_alpha[!is.na(df_alpha$WOMAC_pain),]
  df_alpha[df_alpha=="No response"]<-"NR"
  df_alpha[df_alpha=="Response"]<-"R"
  
  # drop unpaired HostSubjectId
  df_alpha <- df_alpha |>
    filter(n() > 1, .by = HostSubjectId)
  
  # capitalize Timepoint
  df_alpha$Timepoint = str_to_title(df_alpha$Timepoint)
  
  ### Alpha Diversity Statistics ###
  
  # create tables for storing wilcoxon and ttest results
  stats.table.all <- matrix(data = NA, nrow = 1, ncol = 3)
  colnames(stats.table.all) <- c("alpha div", "wilcoxon", "ttest")
  
  # calculate adiv
  stats.table.all[1,1] <- 'shannon_entropy' # colnames(df_alpha)[1] # Timepoint is the first column
  stats.table.all[1,2] <- wilcox.test(shannon_entropy ~ Timepoint, data = df_alpha, paired = TRUE)$p.value
  stats.table.all[1,3] <- t.test(shannon_entropy ~ Timepoint, data = df_alpha, paired = TRUE)$p.value
  
  # save
  ft.all = paste0(outdir,'/',basename(dir),"_paired_alpha_stats.csv")
  write.csv(file = ft.all, stats.table.all)
  
  ### Alpha Diversity Boxplots ###
  
  # variable of interest
  a <- 'shannon_entropy'
  # create filenames
  filename_table = paste(a, basename(dir), "all_table.csv", sep = "_")
  filename_box.plot = paste(a, basename(dir), "all_box.plot.pdf", sep = "_")  
  filename_line.plot = paste(a, basename(dir), "all_line.plot.pdf", sep = "_")  
  
  # create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
  # subset and spread dataset into Timepoint columns
  # then calculate delta relative abundance between post and pre timepoints
  d.div <- df_alpha %>%
    subset(select = c("HostSubjectId", "Timepoint", a)) %>%
    spread(key = "Timepoint", value = a) %>%
    mutate(diff_time.pair = (get('Post') - get('Pre'))) %>%
    mutate(
      Change.type_time.pair = case_when(
        sign(diff_time.pair) == 1 ~ "1_up",
        sign(diff_time.pair) == -1 ~ "2_down",
        TRUE ~ "3_no.change"
      )
    ) 
  
  # merge with original dataset
  d.final <- gather(d.div, 'Pre', 'Post', key = "Timepoint", value = "shannon_entropy")
  
  # rewrite order of factors
  d.final$Timepoint <- factor(d.final$Timepoint, levels = c('Pre', 'Post'))
  
  # save dataset
  # ft = paste(outdir,filename_table, sep = "")
  # write.csv(d.final, file = ft)
  
  # plot boxplot
  p <- ggplot(data = df_alpha, aes(x = Timepoint, y = shannon_entropy, fill = Timepoint)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                 color = "black", size = 0.5, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.5) +
    scale_x_discrete(labels = c('Pre', 'Post')) +
    scale_fill_manual(values = col1) +      
    xlab(NULL) +
    ylab("Shannon Entropy") +
    abkg
  
  # don't plot the boxplot...
  # fpb = paste(outdir, filename_box.plot, sep = "")
  # pdf(file = fpb, height = 4.5, width = 4)
  # plot(p)
  # dev.off()
  
  # plot lineplot
  pv = as.numeric(stats.table.all[1,2])
  pa <- ggplot(data = d.final, aes(x = Timepoint, y = shannon_entropy, fill = Timepoint)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                 color = "black", size = 0.5, width = 0.5) +
    geom_point(shape = 20, size = 3, color = "black") +
    geom_line(data = subset(d.final), aes(group = HostSubjectId, color = Change.type_time.pair, linetype = Change.type_time.pair), size = 0.5) +
    geom_text_repel(data = subset(d.final, Timepoint == 'Pre'), aes(label = HostSubjectId), 
                    nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
    scale_x_discrete(labels = c('Pre', 'Post')) +
    scale_fill_manual(values = col1) +
    scale_linetype_manual(values = line1, name = "Difference", labels = c("Increase", "Decrease", "No change")) +
    scale_color_manual(values = col2, name = "Difference", labels = c("Increase", "Decrease", "No change")) +      
    xlab(NULL) +
    ylab("Shannon Entropy") +
    ggtitle(paste0('p=',as.character(round(pv,3)))) + 
    abkg
  
  if (pv < 0.05){
    pa <- pa + geom_signif(comparisons = list(c("Pre", "Post")), 
                         map_signif_level = TRUE,
                         textsize = 12, # 4,
                         step_increase = 0.05,
                         test = "t.test",
                         vjust = 0.5, # -0.5,
                         annotations = add_significance(pv))
  }
  
  fpl = paste(outdir, filename_line.plot, sep = "")
  pdf(file = fpl, height = 3, width = 4)
  plot(pa)
  dev.off()
  
  ## beta plots
  col1 <- c("#929aab", "#ce2525")
  col1unpaired <- c("green4", "lightblue")
  col2 <- c("#f3a333", "#0074e4", "#8f8787")
  
  
      #theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))
  
  # function to specify that axis labels have 2 decimal places
  f.dec <- function(x){
    format(round(x, 2), nsmall = 2)
  }
  
  # directory for storing files
  outdir = "/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs32/"
  
  # list of distance methods
  dists <- c('bray_curtis')
  
  
  dir=paired_dirs[i]
  # load data
  df <- read.delim(file=paste0(dir,"pcoa.tsv"),
                   row.names=1)
  
  # drop unpaired HostSubjectId
  df <- df |>
    filter(n() > 1, .by = HostSubjectId)
  
  # capitalize Timepoint
  df$Timepoint = str_to_title(df$Timepoint)
  
  # order factors for legend
  df$Timepoint <- factor(df$Timepoint, levels=c('Pre', 'Post'))
  
  for (j in seq_along(dists)) {
    # create filenames
    filename_plot = paste(basename(dir),"bdiv", dists[j], "plot.pdf", sep = "_")
    
    # plot beta diversity
    pv = paired_p[i]
    pb <- ggplot() + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
      geom_point(data = df, aes(x = PC1, y = PC2, color = Timepoint), size=4) + #shape = sib_02),size=4) +
      scale_color_manual(values = c("Pre" = col1[1], "Post" = col1[2])) + #col1, labels = c("Unaffected", "RA")) +
      guides(shape = "none") + 
      ggtitle(paste0('p=',pv)) +
      scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
      scale_y_continuous(labels = f.dec) +
      bbkg # 2 decimal places on y-axis
    
    # save plot
    fp = paste(outdir, filename_plot, sep = "")
    pdf(file = fp, height = 3, width = 4)
    plot(pb)
    dev.off()
  }
  ## barplots
  dir = paired_dirs[i]
  col1 = paired_cols[[i]]
  
  # read in df
  df = read.table(paste0(dir,'barplot_data.tsv'), 
                  sep = '\t', header = TRUE)#, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  # use reorder to reorder the factor
  df$feature = with(df, reorder(feature, logp))
  
  # reorder 
  df$direction = factor(df$direction,c('pre','post'))
  
  # apply the same order to the whole data
  df$feature = factor(df$feature, levels = levels(df$feature))
  # order Enrichment_Significance
  # df$Enrichment_Significance = factor(df$Enrichment_Significance, levels = c('pre_fdr','pre_nofdr','pre_ns','post_ns','post_nofdr','post_fdr','nonresp_fdr','nonresp_nofdr','nonresp_ns','resp_ns','resp_nofdr','resp_fdr'))
  # df$Enrichment_Significance = factor(df$Enrichment_Significance, levels = c('post_fdr','post_nofdr','post_ns','pre_ns','pre_nofdr','pre_fdr','resp_fdr','resp_nofdr','resp_ns','nonresp_ns','nonresp_nofdr','nonresp_fdr'))
  df$Enrichment_Significance = factor(df$Enrichment_Significance, levels = c('Post_FDR','Post_NonFDR','Post_NS','Pre_NS','Pre_NonFDR','Pre_FDR',
                                                                             'R_FDR','R_NonFDR','R_NS','NR_NS','NR_NonFDR','NR_FDR'))
  
  pc <- ggplot(df, aes(x = logp, y = feature, fill = Enrichment_Significance)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "sgn(enrichment)*log10p", y = "Feature") +
    scale_fill_manual(values = col1) +
    
    theme_minimal()
  
  fbp = paste0(dir, 'es_barplot.pdf')
  pdf(file = fbp, height = 5, width = 6)
  plot(pc)
  dev.off()
  
  # figure <- ggarrange(pa, pc, pb,
  #                     labels = c("i", "iii", "ii"),
  #                     ncol = 2, nrow = 2)
  # library(cowplot)
  # # figure <- plot_grid(pa, pc, pb, labels = c('i', 'iii', 'ii'), ncol=2,nrow=2)
  # figure <- plot_grid(pa, pb, pc, labels = c('i', 'ii', 'iii'), 
  #                     ncol=3,nrow=1,
  #                     rel_heights = c(1,1,1),
  #                     rel_widths = c(1.25,1.5,2)) + ggtitle(basename(dir))
  # fbp = paste0(dir, 'fullplot.pdf')
  # pdf(file = fbp, height = 4, width = 16)
  # plot(figure)
  # dev.off()
  
  test[[3*(i-1) + 1]] = pa
  test[[3*(i-1) + 2]] = pb
  test[[3*(i-1) + 3]] = pc
}


for (k in 1:3){
  figure <- plot_grid(test[[3*(k-1)+1]], test[[3*(k-1)+2]], test[[3*(k-1)+3]],
                      test[[3*(k-1)+19]], test[[3*(k-1)+20]], test[[3*(k-1)+21]],
                      test[[3*(k-1)+37]], test[[3*(k-1)+38]], test[[3*(k-1)+39]],
                      test[[3*(k-1)+10]], test[[3*(k-1)+11]], test[[3*(k-1)+12]],
                      test[[3*(k-1)+28]], test[[3*(k-1)+29]], test[[3*(k-1)+30]],
      
                      labels = c('A', '', '',
                                  'B', '', '',
                                  'C', '', '',
                                  'D', '', '',
                                  'E', '', ''), 
                      ncol=3,nrow=5,
                      rel_heights = rep(c(1,1,1),5),
                      rel_widths = rep(c(1.25,1.5,2),5)) + ggtitle(basename(paired_dirs[k]))
  fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/figure_',as.character(k),'_fullplot.pdf')
  pdf(file = fbp, height = 20, width = 16)
  plot(figure)
  dev.off()
}




# 15 iterations
# 3 figures come out of this
# each one contains 5 panels each with 3 subpanels = 45 total panels
# each loop creates one panel 
# iterations i, i + 3, i + 6, i + 9 and i + 12 belong to the same figure
# to store this without overlapping, we create a list of 45 
# l[0+i], l[15+i], l[30 + i]
# when we construct the figure, we will iterate via 





test = list()
### unpaired
for (i in 1:length(unpaired_dirs)){
  # alpha plots
  # background theme
  
    # choose colors
  col1 <- c("#929aab", "#ce2525")
  # col1unpaired <- c("green4", "purple")
  col1unpaired <- c("green4", "lightblue")
  col2 <- c("#f3a333", "#0074e4", "#8f8787")
  
  # choose line types
  line1 <- c("solid", "dashed", "dotted")
  
  # set output directory
  outdir="~/Desktop/clemente_lab/Projects/oa/outputs/jobs31/"
  
  #for (i in 2){
  dir = unpaired_dirs[i]
  print(dir)
  
  df_alpha = read.table(paste0(dir,'metadata.tsv'), 
                        sep = '\t', 
                        header = TRUE, 
                        row.names = 1,
                        check.names = FALSE, 
                        na.strings = "NA")
  
  #drop na
  df_alpha = df_alpha[!is.na(df_alpha$WOMAC_pain),]
  df_alpha[df_alpha=="No response"]<-"NR"
  df_alpha[df_alpha=="Response"]<-"R"
  
  ### Alpha Diversity Statistics ###
  
  # create tables for storing wilcoxon and ttest results
  stats.table.all <- matrix(data = NA, nrow = 1, ncol = 3)
  colnames(stats.table.all) <- c("alpha div", "wilcoxon", "ttest")
  
  # calculate adiv
  stats.table.all[1,1] <- 'shannon_entropy' # colnames(df_alpha)[1] # Timepoint is the first column
  stats.table.all[1,2] <- wilcox.test(shannon_entropy ~ WOMAC_P_Response, data = df_alpha, paired = FALSE)$p.value
  stats.table.all[1,3] <- t.test(shannon_entropy ~ WOMAC_P_Response, data = df_alpha, paired = FALSE)$p.value
  
  # save
  ft.all = paste0(outdir,'/',basename(dir),"_unpaired_alpha_stats.csv")
  write.csv(file = ft.all, stats.table.all)
  
  ### Alpha Diversity Boxplots ###
  
  # variable of interest
  a <- 'shannon_entropy'
  # create filenames
  filename_table = paste(a, basename(dir), "all_table.csv", sep = "_")
  filename_box.plot = paste(a, basename(dir), "all_box.plot.pdf", sep = "_")  
  
  # rewrite order of factors
  df_alpha$WOMAC_P_Response <- factor(df_alpha$WOMAC_P_Response, levels = c('R', 'NR'))
  
  # plot boxplot
  pv = as.numeric(stats.table.all[1,2])
  pa <- ggplot(data = df_alpha, aes(x = WOMAC_P_Response, y = shannon_entropy, fill = WOMAC_P_Response)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                 color = "black", size = 0.5, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.5) +
    scale_x_discrete(labels = c('R', 'NR')) +
    scale_fill_manual(values = col1unpaired) +      
    xlab(NULL) + 
    ggtitle(paste0('p=',as.character(round(pv,3)))) +  # https://stackoverflow.com/questions/39335005/add-p-value-and-r-on-ggplot-follow-up
    ylab("Shannon Entropy") +
    abkg
  
  # add significance
  if (pv < 0.05){
    pa <- pa + geom_signif(comparisons = list(c("NR", "R")), 
                         map_signif_level = TRUE,
                         textsize = 12, # 4,
                         step_increase = 0.05,
                         test = "wilcox.test",
                         vjust = 0.5, # -0.5,
                         annotations = add_significance(pv))
  }
  
  fpb = paste(outdir, filename_box.plot, sep = "")
  pdf(file = fpb, height = 3, width = 4)
  plot(pa)
  dev.off()
  
  
  ## beta plots
  # function to specify that axis labels have 2 decimal places
  f.dec <- function(x){
    format(round(x, 2), nsmall = 2)
  }
  
  # directory for storing files
  outdir = "/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs32/"
  
  # list of distance methods
  dists <- c('bray_curtis')
  
  col1 <- c("#929aab", "#ce2525")
  col1unpaired <- c("green4", "lightblue")
  col2 <- c("#f3a333", "#0074e4", "#8f8787")
  
  dir=unpaired_dirs[i]
  # load data
  df <- read.delim(file=paste0(dir,"pcoa.tsv"),
                   row.names=1)
  
  # dropna
  df = df[!is.na(df$WOMAC_P_Response),]
  df[df=="No response"]<-"NR"
  df[df=="Response"]<-"R"
  
  # order factors for legend
  df$WOMAC_P_Response <- factor(df$WOMAC_P_Response, levels = c('R', 'NR'))
  
  for (j in seq_along(dists)) {
    # create filenames
    filename_plot = paste(basename(dir),"bdiv", dists[j], "plot.pdf", sep = "_")
    
    # plot beta diversity
    pv = unpaired_p[i]
    pb <- ggplot() + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
      geom_point(data = df, aes(x = PC1, y = PC2, color = WOMAC_P_Response), size=4) + #shape = sib_02),size=4) +
      scale_color_manual(values = c("R" = col1unpaired[1], "NR" = col1unpaired[2])) + #col1, labels = c("Unaffected", "RA")) +
      guides(shape = "none") + 
      ggtitle(paste0('p=',pv)) +
      scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
      scale_y_continuous(labels = f.dec) +   # 2 decimal places on y-axis
      bbkg 

    # save plot
    fp = paste(outdir, filename_plot, sep = "")
    pdf(file = fp, height = 3, width = 4)
    plot(pb)
    dev.off()
  }
  
  # barplots
  dir = unpaired_dirs[i]
  col1 = unpaired_cols[[i]]
  
  # read in df
  df = read.table(paste0(dir,'barplot_data.tsv'), 
                  sep = '\t', header = TRUE)#, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  # use reorder to reorder the factor
  df$feature = with(df, reorder(feature, logp))
  
  # reorder 
  df$direction = factor(df$direction,c('NR','R'))
  
  # apply the same order to the whole data
  df$feature = factor(df$feature, levels = levels(df$feature))
  
  # order Enrichment_Significance
  df$Enrichment_Significance = factor(df$Enrichment_Significance, levels = c('Post_FDR','Post_NonFDR','Post_NS','Pre_NS','Pre_NonFDR','Pre_FDR',
                                                                             'R_FDR','R_NonFDR','R_NS','NR_NS','NR_NonFDR','NR_FDR'))
  
  pc <- ggplot(df, aes(x = logp, y = feature, fill = Enrichment_Significance)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "sgn(enrichment)*log10p", y = "Feature") +
    scale_fill_manual(values = col1) +
    
    theme_minimal()
  
  fbp = paste0(dir, 'es_barplot.pdf')
  pdf(file = fbp, height = 5, width = 6)
  plot(pc)
  dev.off()

  test[[3*(i-1) + 1]] = pa
  test[[3*(i-1) + 2]] = pb
  test[[3*(i-1) + 3]] = pc
}

for (k in 1:2){
  figure <- plot_grid(test[[3*(k-1)+7]], test[[3*(k-1)+8]], test[[3*(k-1)+9]],
                      test[[3*(k-1)+25]], test[[3*(k-1)+26]], test[[3*(k-1)+27]],
                      test[[3*(k-1)+19]], test[[3*(k-1)+20]], test[[3*(k-1)+21]],
                      test[[3*(k-1)+1]], test[[3*(k-1)+2]], test[[3*(k-1)+3]],
                      test[[3*(k-1)+13]], test[[3*(k-1)+14]], test[[3*(k-1)+15]],
                      
                      labels = c('A', '', '',
                                 'B', '', '',
                                 'C', '', '',
                                 'D', '', '',
                                 'E', '', ''), 
                      ncol=3,nrow=5,
                      rel_heights = rep(c(1,1,1),5),
                      rel_widths = rep(c(1.25,1.5,2),5)) + ggtitle(basename(paired_dirs[k]))
  fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/figure_up_',as.character(k),'_fullplot.pdf')
  pdf(file = fbp, height = 20, width = 16)
  plot(figure)
  dev.off()
}
  

test = list()
for (i in 1:length(diff_dirs)){
  dir = diff_dirs[i]
  col1 = diff_cols[[i]]
  
  # read in df
  df = read.table(paste0(dir,'barplot_data.tsv'), 
                  sep = '\t', header = TRUE)#, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  # use reorder to reorder the factor
  df$feature = with(df, reorder(feature, logp))
  
  # reorder 
  df$direction = factor(df$direction,c('NR','R'))
  
  # order Enrichment_Significance
  # df$Enrichment_Significance = factor(df$Enrichment_Significance, levels = c('pre_fdr','pre_nofdr','pre_ns','post_ns','post_nofdr','post_fdr','nonresp_fdr','nonresp_nofdr','nonresp_ns','resp_ns','resp_nofdr','resp_fdr'))
  df$Enrichment_Significance = factor(df$Enrichment_Significance, levels = c('Post_FDR','Post_NonFDR','Post_NS','Pre_NS','Pre_NonFDR','Pre_FDR',
                                                                             'R_FDR','R_NonFDR','R_NS','NR_NS','NR_NonFDR','NR_FDR'))
  
  # apply the same order to the whole data
  df$feature = factor(df$feature, levels = levels(df$feature))
  
  p <- ggplot(df, aes(x = logp, y = feature, fill = Enrichment_Significance)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "log10p", y = "Feature") +
    scale_fill_manual(values = col1) +
    
    theme_minimal()
  
  fbp = paste0(dir, 'es_barplot.pdf')
  pdf(file = fbp, height = 5, width = 6)
  plot(p)
  dev.off()

  test[[i]] = p
}

figure <- plot_grid(test[[1]],
                    test[[2]],
                    test[[3]],
                    test[[4]],
                    test[[5]],
                    
                    labels = c('A', 'B', 'C','D','E'),
                    ncol=1,nrow=5,
                    rel_heights = c(1,1,1,1,1),
                    rel_widths = c(1,1,1,1,1)) + ggtitle(basename(paired_dirs[k]))
fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/figure_diffs_',as.character(k),'_fullplot.pdf')
pdf(file = fbp, height = 18, width = 8)
plot(figure)
dev.off()



