################################# 
## R script                    ##
## Project: TwinsRA            ##
## Data: Shotgun metagenomics  ##
## Author: KB                  ##
## Date: 5/17/24               ##
#################################

# libraries ---------------------------------------------------------------
library(plyr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(gridtext)
library(caret)

# Load data ---------------------------------------------------------------
setwd('~/Desktop/clemente_lab/Projects/oa/')

for (y in c('saliva', 'stool')){
  if (y == 'stool'){
    inputs = c('outputs/jobs16/data_processing/summary_df_resample_1.txt')
  } else if (y == 'saliva'){
    inputs = c('outputs/jobs15/data_processing/summary_df_resample_1.txt')
  }
  metadata = read.table(file=paste('inputs/df_quant_',y,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  metadata = read.table(file=paste('inputs/cutie_df_abun_',y,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  metadata = read.table(file=paste('inputs/cutie_df_taxa_',y,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  
  # Compile Cutie results and calculate corr difference ----------------
  res_dir = '~/Desktop/clemente_lab/Projects/oa/'
  # cutie_res = {}
  # grab first input
  fp = inputs[1]
  df = read.table(fp, sep='\t', header=T, 
                  stringsAsFactors=F, check.names=F)
  # filters to keep variables that are in the metadata file; make sure no OTU here
  filt = sapply(df[1:2], function(x) x %in% c(colnames(metadata)))
  
  # keep columns 1, 2 and 4 and keep rows where one of the two entries is in the metadata file list of vars
  # keep columns 1, 2 and 4 and keep rows where BOTH of the two entries is in the metadata file list of vars
  # df_filt = df[apply(filt, 2, sum)==1, c(1:2, 4)]
  df_filt = df[apply(filt, 1, sum)==1, c(1:2, 4)]
  corr_matrix = spread(df_filt, var2, correlations) %>%
    'rownames<-'(.[, 1]) %>% .[, -1]
  # cutie_res[[i]] = corr_matrix
  
  # grab corr matrix
  corr_diff = corr_matrix
  
  # subset specific analyses
  #colnames(corr_diff)[c(17, 23, 33)] = c('IFN-g', 'IL-1b', 'TNF-a')
  
  # Heatmap -----------------------------------------------------------------
  # Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
  #corr_diff <- corr_diff[,colSums(is.na(corr_diff))<nrow(corr_diff)]
  #corr_diff <- corr_diff[rowSums(is.na(corr_diff))<ncol(corr_diff),]
  #corr_diff[is.na(corr_diff)] <- 0
  data_matrix = as.matrix(corr_diff) # %>% na.omit
  # is replacing with 0 valid?
  #dend = as.dendrogram(hclust(dist(data_matrix), 
  #                            method='average'))

  # create annotations
  df_anno_cols = read.table(file=paste('inputs/col_',y,'.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  cols_anno = df_anno_cols$var_type
  
  df_anno_rows = read.table(file=paste('inputs/row_',y,'.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  rows_anno = df_anno_rows$var_type
  
  #Create a custom color scale
  # https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
  library(RColorBrewer)
  col_colors <- brewer.pal(length(unique(cols_anno)),"Set1")
  names(col_colors) <- levels(df_anno_cols$var_type)
  col_scale <- scale_colour_manual(name = "Feature Type",values = col_colors)

  row_colors <- brewer.pal(length(unique(rows_anno)),"Set1")
  names(row_colors) <- levels(df_anno_rows$var_type)
  row_scale <- scale_colour_manual(name = "Feature Type",values = row_colors)
  
  # remove nearzerovar
  col_idx <- nearZeroVar(data_matrix)
  if (identical(col_idx, integer(0))){
    data_matrix <- data_matrix
  } else {
    data_matrix <- data_matrix[,-c(col_idx)]
    cols_anno <- cols_anno[-c(col_idx)]
  }
  # set color annotation
  # https://www.biostars.org/p/317349/
  col_ann <- data.frame(cols_anno)
  colnames(col_ann) <- c('Type')#, 'Type2')
  col_colours <- list('Type' = c('Sig' = myColors[1], 
                             'Nonsig' = myColors[2]))
  #                           'Metagenomic_Pathways' = myColors[8]))
  #  
  colAnn <- HeatmapAnnotation(df = col_ann,
                              which = 'col',
                              col = col_colours,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))
  
  row_ann <- data.frame(rows_anno)
  colnames(row_ann) <- c('Type')#, 'Type2')
  row_colours <- list('Type' = c('BA' = myColors[1], 
                             'NonBA' = myColors[2]))
  #                           'Metagenomic_Pathways' = myColors[8]))
  #  
  rowAnn <- HeatmapAnnotation(df = row_ann,
                              which = 'row',
                              col = row_colours,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))
  
  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs17/heatmap_corr_',y,'.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 8, height=8)  # Image size
  hm <- Heatmap(data_matrix, use_raster = FALSE, name='Spearman',
                width = unit(8, "cm"),
                row_dend_width = unit(1.5, "cm"),
                
                clustering_method_rows = "average",  # UPGMA
                clustering_method_columns = "average",
                
                row_title = 'Metabolites',
                row_names_gp = gpar(fontsize = 4),
                show_row_names = FALSE,
                
                column_title = 'Taxa',
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 4),
                
                top_annotation = colAnn, #columnAnnotation(VarType=cols_anno, col=colScale),
                left_annotation = rowAnn,
                #heatmap_legend_param = list(
                #  title = title,
                #  direction = "horizontal"
                #)
  )
  
  hm_draw = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="right"
  )
  dev.off()
}
