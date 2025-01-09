# Load required library
library(ggplot2)


# specify directories
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



for (i in 1:length(paired_dirs)){
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
  # order combhue
  # df$combhue = factor(df$combhue, levels = c('pre_fdr','pre_nofdr','pre_ns','post_ns','post_nofdr','post_fdr','nonresp_fdr','nonresp_nofdr','nonresp_ns','resp_ns','resp_nofdr','resp_fdr'))
  df$combhue = factor(df$combhue, levels = c('post_fdr','post_nofdr','post_ns','pre_ns','pre_nofdr','pre_fdr','resp_fdr','resp_nofdr','resp_ns','nonresp_ns','nonresp_nofdr','nonresp_fdr'))
  
  p <- ggplot(df, aes(x = logp, y = feature, fill = combhue)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "log10p", y = "Feature") +
    scale_fill_manual(values = col1) +
    
    theme_minimal()
  
  fbp = paste0(dir, 'es_barplot.pdf')
  pdf(file = fbp, height = 5, width = 6)
  plot(p)
  dev.off()
}

for (i in 1:length(unpaired_dirs)){
  dir = unpaired_dirs[i]
  col1 = unpaired_cols[[i]]
  
  # read in df
  df = read.table(paste0(dir,'barplot_data.tsv'), 
                  sep = '\t', header = TRUE)#, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  # use reorder to reorder the factor
  df$feature = with(df, reorder(feature, logp))
  
  # reorder 
  df$direction = factor(df$direction,c('No response','Response'))
  
  # apply the same order to the whole data
  df$feature = factor(df$feature, levels = levels(df$feature))
  
  # order combhue
  # df$combhue = factor(df$combhue, levels = c('pre_fdr','pre_nofdr','pre_ns','post_ns','post_nofdr','post_fdr','nonresp_fdr','nonresp_nofdr','nonresp_ns','resp_ns','resp_nofdr','resp_fdr'))
  df$combhue = factor(df$combhue, levels = c('post_fdr','post_nofdr','post_ns','pre_ns','pre_nofdr','pre_fdr','resp_fdr','resp_nofdr','resp_ns','nonresp_ns','nonresp_nofdr','nonresp_fdr'))
  
  p <- ggplot(df, aes(x = logp, y = feature, fill = combhue)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "log10p", y = "Feature") +
    scale_fill_manual(values = col1) +
    
    theme_minimal()
  
  fbp = paste0(dir, 'es_barplot.pdf')
  pdf(file = fbp, height = 5, width = 6)
  plot(p)
  dev.off()
}


for (i in 1:length(diff_dirs)){
  dir = diff_dirs[i]
  col1 = diff_cols[[i]]
  
  # read in df
  df = read.table(paste0(dir,'barplot_data.tsv'), 
                  sep = '\t', header = TRUE)#, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  # use reorder to reorder the factor
  df$feature = with(df, reorder(feature, logp))
  
  # reorder 
  df$direction = factor(df$direction,c('No response','Response'))
  
  # order combhue
  # df$combhue = factor(df$combhue, levels = c('pre_fdr','pre_nofdr','pre_ns','post_ns','post_nofdr','post_fdr','nonresp_fdr','nonresp_nofdr','nonresp_ns','resp_ns','resp_nofdr','resp_fdr'))
  df$combhue = factor(df$combhue, levels = c('post_fdr','post_nofdr','post_ns','pre_ns','pre_nofdr','pre_fdr','resp_fdr','resp_nofdr','resp_ns','nonresp_ns','nonresp_nofdr','nonresp_fdr'))
  
  # apply the same order to the whole data
  df$feature = factor(df$feature, levels = levels(df$feature))
  
  p <- ggplot(df, aes(x = logp, y = feature, fill = combhue)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "log10p", y = "Feature") +
    scale_fill_manual(values = col1) +
    
    theme_minimal()
  
  fbp = paste0(dir, 'es_barplot.pdf')
  pdf(file = fbp, height = 5, width = 6)
  plot(p)
  dev.off()
}
