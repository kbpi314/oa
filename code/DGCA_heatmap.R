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
  metadata = read.table(file=paste('inputs/cutie_df_abun_',y,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  
  # Compile Cutie results and calculate corr difference ----------------
  res_dir = '~/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19'
  # cutie_res = {}
  # grab first input
  fp = inputs[1]
  df = read.table(fp, sep='\t', header=T, 
                  stringsAsFactors=F, check.names=F)
  # filters to keep variables that are in the metadata file; make sure no OTU here
  filt = sapply(df[1:2], function(x) x %in% c(colnames(metadata)))
  
  # keep columns 1, 2 and 4 and keep rows where one of the two entries is in the metadata file list of vars
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
  corr_diff <- corr_diff[rowSums(is.na(corr_diff))<ncol(corr_diff),]
  corr_diff[is.na(corr_diff)] <- 0
  data_matrix = as.matrix(corr_diff) # %>% na.omit
  # is replacing with 0 valid?
  #dend = as.dendrogram(hclust(dist(data_matrix), 
  #                            method='average'))

  # create annotations
  #df_anno = read.table(file=paste('inputs/df_quant_R_corr_labels.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  #cols_anno = df_anno$var_type
  
  #Create a custom color scale
  # https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
  library(RColorBrewer)
  myColors <- brewer.pal(8,"Set1")
  #names(myColors) <- levels(df_anno$var_type)
  colScale <- scale_colour_manual(name = "Feature Type",values = myColors)
  
  # remove nearzerovar
  col_idx <- nearZeroVar(data_matrix)
  if (identical(col_idx, integer(0))){
    data_matrix <- data_matrix
  } else {
    data_matrix <- data_matrix[,-c(col_idx)]
    # cols_anno <- cols_anno[-c(col_idx)]
  }
  # set color annotation
  # https://www.biostars.org/p/317349/
  #ann <- data.frame(cols_anno)
  #colnames(ann) <- c('Type')#, 'Type2')
  #colours <- list('Type' = c('Clinical' = myColors[1], 
  #                           'Plasma_Autoantibodies' = myColors[2],
  #                           'Fecal_Autoantibodies' = myColors[3],
  #                           'Plasma_Cytokines' = myColors[4],
  #                           'Serum_Fatty_Acids' = myColors[5],
  #                           'Stool_Fatty_Acids' = myColors[6],
  #                           'Metagenomic_Taxa' = myColors[7],
  #                           'Metagenomic_Pathways' = myColors[8]))
  #  
  #colAnn <- HeatmapAnnotation(df = ann,
  ##                            which = 'col',
  #                            col = colours,
  #                            annotation_width = unit(c(1, 4), 'cm'),
  #                            gap = unit(1, 'mm'))
  
  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs17/heatmap_corr_',y,'.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 8, height=8)  # Image size
  hm <- Heatmap(data_matrix, use_raster = FALSE,
                width = unit(8, "cm"),
                row_dend_width = unit(1.5, "cm"),
                
                clustering_method_rows = "average",  # UPGMA
                clustering_method_columns = "average",
                
                # row_title = y,
                row_names_gp = gpar(fontsize = 4),
                show_row_names = FALSE,
                
                # column_title = 'Quant Vars',
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 4),
                
                # top_annotation = colAnn, #columnAnnotation(VarType=cols_anno, col=colScale),
                
                #heatmap_legend_param = list(
                #  title = title,
                #  direction = "horizontal"
                #)
  )
  
  hm = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="right"
  )
  dev.off()
}
