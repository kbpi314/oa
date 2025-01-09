################################# 
## R script                    ##
## Project: Twinspost            ##
## Alpha diversity             ##
## Data: Shotgun metagenomics  ##
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

############################################################################
############################################################################
############################################################################

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 24, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 18))+#, color = "black")) +
  theme(axis.title.y = element_text(size = 24, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black"))

# choose colors
col1 <- c("#929aab", "#ce2525")
# col1unpaired <- c("green4", "purple")
col1unpaired <- c("green4", "lightblue")
col2 <- c("#f3a333", "#0074e4", "#8f8787")

# choose line types
line1 <- c("solid", "dashed", "dotted")

# set output directory
outdir="~/Desktop/clemente_lab/Projects/oa/outputs/jobs31/"

for (i in 1:length(paired_dirs)){
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
  
  # drop unpaired HostSubjectId
  df_alpha <- df_alpha |>
      filter(n() > 1, .by = HostSubjectId)
  
  # capitalize Timepoint
  df_alpha$Timepoint = toTitleCase(df_alpha$Timepoint)
  
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
    bkg
  
  # don't plot the boxplot...
  # fpb = paste(outdir, filename_box.plot, sep = "")
  # pdf(file = fpb, height = 4.5, width = 5)
  # plot(p)
  # dev.off()
  
  # plot lineplot
  pv = as.numeric(stats.table.all[1,2])
  p <- ggplot(data = d.final, aes(x = Timepoint, y = shannon_entropy, fill = Timepoint)) +
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
    bkg + 
    theme(title = element_text(size = 18, color = "black")) 
  
  if (pv < 0.05){
    p <- p + geom_signif(comparisons = list(c("Pre", "Post")), 
                         map_signif_level = TRUE,
                         textsize = 12, # 4,
                         step_increase = 0.05,
                         test = "t.test",
                         vjust = 0.5, # -0.5,
                         annotations = add_significance(pv))
  }
  
  fpl = paste(outdir, filename_line.plot, sep = "")
  pdf(file = fpl, height = 6, width = 8)
  plot(p)
  dev.off()
}

bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 18, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 18))+#, color = "black")) +
  theme(axis.title.y = element_text(size = 24, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  theme(legend.title = element_blank())
  # theme(legend.title = element_text(size = 18, face = "bold", color = "black"))

### unpaired
for (i in 1:length(unpaired_dirs)){
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
  df_alpha$WOMAC_P_Response <- factor(df_alpha$WOMAC_P_Response, levels = c('Response', 'No response'))
  
  # plot boxplot
  pv = as.numeric(stats.table.all[1,2])
  p <- ggplot(data = df_alpha, aes(x = WOMAC_P_Response, y = shannon_entropy, fill = WOMAC_P_Response)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                 color = "black", size = 0.5, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.5) +
    scale_x_discrete(labels = c('Response', 'No response')) +
    scale_fill_manual(values = col1unpaired) +      
    xlab(NULL) + 
    ggtitle(paste0('p=',as.character(round(pv,3)))) +  # https://stackoverflow.com/questions/39335005/add-p-value-and-r-on-ggplot-follow-up
    ylab("Shannon Entropy") +
    bkg + 
    theme(title = element_text(size = 18, color = "black")) 
  
  # add significance
  if (pv < 0.05){
    p <- p + geom_signif(comparisons = list(c("No response", "Response")), 
                         map_signif_level = TRUE,
                         textsize = 12, # 4,
                         step_increase = 0.05,
                         test = "wilcox.test",
                         vjust = 0.5, # -0.5,
                         annotations = add_significance(pv))
  }
  
  fpb = paste(outdir, filename_box.plot, sep = "")
  pdf(file = fpb, height = 4.5, width = 8)
  plot(p)
  dev.off()
}


