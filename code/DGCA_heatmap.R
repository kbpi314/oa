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

for (y in c('stool', 'saliva')){ # saliva
  # R and then NR
  if (y == 'stool'){
    inputs = c('outputs/jobs22/data_processing/summary_df_resample_1.txt', 'outputs/jobs24/data_processing/summary_df_resample_1.txt')
  }
  else if (y == 'saliva'){
    inputs = c('outputs/jobs21/data_processing/summary_df_resample_1.txt', 'outputs/jobs23/data_processing/summary_df_resample_1.txt')
  }
  
  # metadata = read.table(file=paste('inputs/df_quant_',y,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  # metadata = read.table(file=paste('inputs/cutie_df_abun_',y,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
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
  
  # is replacing with 0 valid?
  #dend = as.dendrogram(hclust(dist(data_matrix), 
  #                            method='average'))

  # create annotations
  df_anno_cols = read.table(file=paste('inputs/col_',y,'.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  cols_anno = df_anno_cols$var_type
  
  df_anno_rows = read.table(file=paste('inputs/row_',y,'.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  rows_anno = df_anno_rows$var_type
  
  # split corr_diff matrix into bile acids and non bile acids
  BA_corr_diff <- corr_diff %>% filter(row.names(corr_diff) %in% df_anno_rows[df_anno_rows$var_type == 'BA',]$index) #%>% na.omit
  BA_corr_diff <- BA_corr_diff[,colSums(is.na(BA_corr_diff))<nrow(BA_corr_diff)]
  BA_matrix <- as.matrix(BA_corr_diff)

  non_BA_corr_diff <- corr_diff %>% filter(row.names(corr_diff) %in% df_anno_rows[df_anno_rows$var_type == 'NonBA',]$index) #%>% na.omit
  non_BA_corr_diff <- non_BA_corr_diff[,colSums(is.na(non_BA_corr_diff))<nrow(non_BA_corr_diff)]
  non_BA_matrix <- as.matrix(non_BA_corr_diff)
  
  
  # repeat for nr
  fp = inputs[2]
  df = read.table(fp, sep='\t', header=T, 
                  stringsAsFactors=F, check.names=F)

  filt = sapply(df[1:2], function(x) x %in% c(colnames(metadata)))
  
  df_filt = df[apply(filt, 1, sum)==1, c(1:2, 4)]
  corr_matrix = spread(df_filt, var2, correlations) %>%
    'rownames<-'(.[, 1]) %>% .[, -1]
  
  # grab corr matrix
  corr_diff = corr_matrix
  
  # split corr_diff matrix into bile acids and non bile acids
  BA_corr_diff <- corr_diff %>% filter(row.names(corr_diff) %in% df_anno_rows[df_anno_rows$var_type == 'BA',]$index) #%>% na.omit
  BA_corr_diff <- BA_corr_diff[,colSums(is.na(BA_corr_diff))<nrow(BA_corr_diff)]
  BA_matrix2 <- as.matrix(BA_corr_diff)
  
  non_BA_corr_diff <- corr_diff %>% filter(row.names(corr_diff) %in% df_anno_rows[df_anno_rows$var_type == 'NonBA',]$index) #%>% na.omit
  non_BA_corr_diff <- non_BA_corr_diff[,colSums(is.na(non_BA_corr_diff))<nrow(non_BA_corr_diff)]
  non_BA_matrix2 <- as.matrix(non_BA_corr_diff)
  
  
  
  
  
  # data_matrix = as.matrix(corr_diff) # 
  
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
  # col_idx <- nearZeroVar(data_matrix)
  # if (identical(col_idx, integer(0))){
  #   data_matrix <- data_matrix
  # } else {
  #   data_matrix <- data_matrix[,-c(col_idx)]
  #   cols_anno <- cols_anno[-c(col_idx)]
  # }
  # set color annotation
  # https://www.biostars.org/p/317349/
  col_ann <- data.frame(cols_anno)
  colnames(col_ann) <- c('Type')#, 'Type2')
  col_colours <- list('Type' = c('Sig' = col_colors[1], 
                             'Nonsig' = col_colors[2]))
  #                           'Metagenomic_Pathways' = myColors[8]))
  #  
  colAnn <- HeatmapAnnotation(df = col_ann,
                              which = 'col',
                              col = col_colours,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))
  
  #row_ann <- data.frame(rows_anno)
  #colnames(row_ann) <- c('Type')#, 'Type2')
  #row_colours <- list('Type' = c('BA' = row_colors[1], 
  #                           'NonBA' = row_colors[2]))
  #                           'Metagenomic_Pathways' = myColors[8]))
  #  
  # rowAnn <- HeatmapAnnotation(df = row_ann,
  #                             which = 'row',
  #                            col = row_colours,
  #                            annotation_width = unit(c(1, 4), 'cm'),
  #                            gap = unit(1, 'mm'))
  
  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs17/heatmap_corr_',y,'_r.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 8, height=8)  # Image size
  hm_BA <- Heatmap(BA_matrix, use_raster = FALSE, name='Spearman',
                width = unit(8, "cm"),
                row_dend_width = unit(1.5, "cm"),
                
                clustering_method_rows = "average",  # UPGMA
                clustering_method_columns = "average",
                
                row_title = 'Bile Acids',
                row_names_gp = gpar(fontsize = 4),
                show_row_names = FALSE,
                
                column_title = 'Taxa',
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 4),
                
                # top_annotation = colAnn, #columnAnnotation(VarType=cols_anno, col=colScale),
                # left_annotation = rowAnn,
                #heatmap_legend_param = list(
                #  title = title,
                #  direction = "horizontal"
                #)
  )
  hm_nBA <- Heatmap(non_BA_matrix, use_raster = FALSE, name='Spearman_Non_Bile_Acids',
                   width = unit(8, "cm"),
                   row_dend_width = unit(1.5, "cm"),
                   
                   clustering_method_rows = "average",  # UPGMA
                   clustering_method_columns = "average",
                   
                   row_title = 'Non-Bile Acids',
                   row_names_gp = gpar(fontsize = 4),
                   show_row_names = FALSE,
                   
                   column_title = 'Taxa',
                   show_column_names = FALSE,
                   column_names_gp = gpar(fontsize = 4),
                   
                   #top_annotation = colAnn, #columnAnnotation(VarType=cols_anno, col=colScale),
                   # left_annotation = rowAnn,
                   #heatmap_legend_param = list(
                   #  title = title,
                   #  direction = "horizontal"
                   #)
                   show_heatmap_legend = FALSE
  )
  hm = hm_BA %v% hm_nBA 
  
  hm = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="right"
  )
  dev.off()

  # Draw other set of heatmaps  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs17/heatmap_corr_',y,'_nr.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 8, height=8)  # Image size
  
  hm_BA <- Heatmap(BA_matrix2, use_raster = FALSE, name='Spearman',
                   width = unit(8, "cm"),
                   row_dend_width = unit(1.5, "cm"),
                   
                   row_title = 'Bile Acids',
                   row_names_gp = gpar(fontsize = 4),
                   show_row_names = FALSE,
                   
                   column_title = 'Taxa',
                   show_column_names = FALSE,
                   column_names_gp = gpar(fontsize = 4),
                   
                   row_order = row_order(hm)$Spearman,
                   # column_order = column_order(hm)
                   
  )
  hm_nBA <- Heatmap(non_BA_matrix2, use_raster = FALSE, name='Spearman_Non_Bile_Acids',
                    width = unit(8, "cm"),
                    row_dend_width = unit(1.5, "cm"),

                    row_title = 'Non-Bile Acids',
                    row_names_gp = gpar(fontsize = 4),
                    show_row_names = FALSE,
                    
                    column_title = 'Taxa',
                    show_column_names = FALSE,
                    column_names_gp = gpar(fontsize = 4),
                    
                    row_order = row_order(hm)$Spearman_Non_Bile_Acids,
                    # column_order = column_order(hm),
                    
                    show_heatmap_legend = FALSE
  )
  hm2 = hm_BA %v% hm_nBA 
  
  hm2 = draw(
    hm2, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="right"
  )

  
  
  dev.off()
}

