# import libraries
library(dplyr); library(plyr); library(glue)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)  # colorRamp2()
library(GetoptLong)
library(ggsci)  # pal_jco()
library(gridtext)
options(stringsAsFactors=F)
library(IRdisplay)  # display() for pretty display of R outputs
options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)  # View more table cells (Default: 60, 20)


# load data
fp = "cutie_df.txt"
cutie_df = read.table(fp, sep='\t', header=T, stringsAsFactors=F,
                      check.names=F, quote='', comment.char='')
fp = "lefse_df.txt"
lefse_df = read.table(fp, sep='\t', header=T, check.names=F,
                      quote='', comment.char='')



# Reformat Metaphlan full taxonomy strings
rename_taxa = function(x) {
  strsplit(as.character(x), ';')[[1]] %>% 
    .[length(.)] %>% 
    {gsub('^s__', '', .)} %>%
    setNames(x)
}

# Will convert Cutie output from long to wide format
# p, r, cutie, var1, var2: Column names in Cutie df,
#    replace with your Cutie columns of interest
create_cutie_matrices = function(cutie.df) {
  output = lapply(list(p='p', r='r', c='cutie'), function(v) {
    cutie.df[, c('var1', 'var2', v)] %>%
      reshape(idvar='var1', timevar='var2', direction='wide') %>%
      'rownames<-'(.$var1) %>% .[,-1] %>%
      'colnames<-'(gsub(glue("^{v}[.]"), "", colnames(.))) %>%
      .[order(rownames(.)), order(colnames(.))] %>%
      as.matrix(check.names=F)
  })
  return(output)
}

# Converts output from create_cutie_matrices into ComplexHeatmap function input
create_heatmap_matrices = function(cutie.matrices, p.filt=T, 
                                   lefse.filt=NULL, r.filt=NULL) {
  # Filter function
  filter_df = function(key, df, rowfilt) {
    if (sum(rowfilt)>0) {
      df = df[rowfilt,]
    } else {
      print(glue("{key} filter error: Zero row features"))
    }
    return(df)
  }
  # Filter by correlation strength
  if (!is.null(r.filt)) {
    rowfilt = apply(cutie.matrices$r, 1, function(x) any(abs(x)>=r.filt, na.rm=T))
    print(glue("r filt: {sum(rowfilt)}/{length(rowfilt)} features"))
    cutie.matrices = lapply(cutie.matrices, function(x) filter_df('r', x, rowfilt))
  }
  # Filter by significant Lefse features
  if (!is.null(lefse.filt)) {
    rowfilt = rownames(cutie.matrices$r) %in% lefse.filt
    print(glue("Lefse filt: {sum(rowfilt)}/{length(rowfilt)} features"))
    cutie.matrices = lapply(cutie.matrices, function(x) filter_df('lefse', x, rowfilt))
  }
  # Filter significant correlations
  check.rowfilt = c()
  if (p.filt) {
    rowfilt = apply(cutie.matrices$p, 1, function(x) any(x<0.05, na.rm=T))
    cutie.matrices = lapply(cutie.matrices, function(x) filter_df('pval', x, rowfilt))
  }
  # Omit all NA rows/cols
  rowfilt = apply(cutie.matrices$r, 1, function(x) all(is.na(x)))
  print(glue("NA filt: {sum(rowfilt)}/{length(rowfilt)} features"))
  colfilt = apply(cutie.matrices$r, 2, function(x) !all(is.na(x)))
  if (sum(rowfilt)>0 & all(lapply(colfilt, sum)>0)) {
    output = lapply(cutie.matrices, function(x) x[rowfilt, colfilt[[i]]])
  } else {
    print(glue("Omit NAs filter error: (row={sum(rowfilt)}, col={sum(colfilt)}) features"))
    output = cutie.matrices
  }
  print(glue("Total: {nrow(cutie.matrices$r)} -> {nrow(output$r)} features"))      
  
  return(output)
}


# Plot heatmap
plot_heatmap = function(hm.inputs, lefse.df,
                        fp, f.w, f.h, plot.title=NULL,
                        hm.title, hm.w, hm.h, 
                        annot.r.lefse.ymax,
                        hm.rownames=NULL, hm.colnames=NULL,
                        pch.size = 6,
                        row.km=1, col.km=1, km.gap=1.5) {
  # Annotations
  annot_top = HeatmapAnnotation(
    `Top Annotation` = annotmap_top[colnames(hm_mat$r)],
    col = list(`Top Annotation`=annotcol_top),
    simple_anno_size = unit(3, "mm"), 
    show_annotation_name = F
  )
  annotmap_right = lefse.df[match(rownames(hm_mat$r), lefse.df$feature), ]
  annot_right = rowAnnotation(
    `Right Annotation` = anno_simple(
      annotmap_right$class,
      col = annotcol_right,
      width = unit(3,'mm')
    ),
    Lefse = anno_barplot(
      round(annotmap_right$log_LDA, digits=1),
      gp=gpar(fill=annotcol_right[annotmap_right$class], col=NA),
      width=unit(1.5, 'cm'),
      bar_width=0.8, 
      add_numbers=TRUE,
      axis_param = list(labels_rot=0),
      numbers_rot = 0,
      numbers_gp = gpar(col='white', cex=0.8),
      numbers_offset = unit(-0.55, 'cm'),
      ylim = c(0, annot.r.lefse.ymax)  # transposed
    ),
    annotation_label = c('', gt_render("log10(LDA)", gp = gpar(cex=0.9))),
    gap = unit(1.7, "mm")  # Between 2 annotations
  )
  # Legends
  pch.shape.p = 0
  pch.shape.c = 20
  pch.color = 'black'
  pch.lwd = 1.7
  lgd1 = Legend(  # anno_simple() legend has to be manually created
    title = 'Right Annotation',
    labels = c('Psychotic Case', 'Healthy Control'),
    legend_gp = gpar(fill = annotcol_right)
  )
  lgd2 = Legend(
    title = 'LEfSe Annotation',
    labels = c('Psychotic Case', 'Healthy Control'),
    legend_gp = gpar(fill = annotcol_right)
  )
  hm_color = colorRamp2(
    c(-1,-0.5,0,0.5,1), 
    c("turquoise4", 'turquoise3', "white", "sienna1", 'sienna2')
  )  # Alternative color scheme: c("dodgerblue3", 'dodgerblue1', "white", "brown1", 'brown3')
  lgd3 = Legend(
    title="Correlation r",
    direction='horizontal',
    col_fun=hm_color,
    at=c(-1,0,1), 
    legend_width=unit(2.5, "cm"),
    grid_height=unit(3, "mm")
  )
  lgd4 = Legend(
    title = "Significant correlations",
    labels = c('P < 0.05', 'True Positive'),
    type = "points", 
    pch = c(pch.shape.p, pch.shape.c),
    legend_gp = gpar(lwd=1.5, col=1), 
    background = NULL, size = unit(c(4,3), "mm")
  )
  lgd = packLegend(lgd1, lgd2, lgd3, lgd4)
  # Heatmaps
  ht_opt(
    heatmap_row_names_gp = gpar(fontsize=11),
    heatmap_column_names_gp = gpar(fontsize=11),
    heatmap_column_title_gp = gpar(fontsize=11),
    ROW_ANNO_PADDING=unit(2.7, 'mm'),
    COLUMN_ANNO_PADDING=unit(2, 'mm')
  )
  hm.rownames = if (is.null(hm.rownames)) rownames(hm.inputs$r) else hm.rownames
  hm.colnames = if (is.null(hm.colnames)) colnames(hm.inputs$r) else hm.colnames
  hm_main = Heatmap(
    hm.inputs$r, col = hm_color, na_col = 'black',
    column_title = hm.title,        
    column_title_gp = gpar(fontsize = 12),
    column_km = col.km,
    column_gap = unit(km.gap, "mm"),
    row_labels = gt_render(
      hm.rownames, padding = unit(c(0,0,0,1.2), "mm") # top,right,bottom,left
    ),
    column_labels = gt_render(hm.colnames, padding=unit(1.2, 'mm')),
    right_annotation = annot_right,
    top_annotation = annot_top,
    width = unit(hm.w, "cm"),
    show_heatmap_legend = F,
    layer_fun = function(j, i, x, y, width, height, fill) {
      v = pindex(hm.inputs$p, i, j)
      l = v<0.05
      if (any(l, na.rm=T)) {
        grid.points(x[l], y[l], pch=pch.shape.p, size=unit(pch.size, "mm"),
                    gp=gpar(lwd=pch.lwd, col=pch.color))
      }
      v = pindex(hm.inputs$c, i, j)
      l = v==1
      if (any(l, na.rm=T)) {
        grid.points(x[l], y[l], pch=pch.shape.c, size=unit(pch.size-3.5, "mm"),
                    gp=gpar(lwd=pch.lwd, col=pch.color))
      }
    }
  )
  # Full heatmap
  png(qq(fp), width=f.w, height=f.h, units='in', res=300)
  draw(
    hm_main, 
    height = unit(hm.h, "cm"),
    ht_gap = unit(2.5, "mm"),
    row_km = row.km, row_gap = unit(km.gap, "mm"),
    column_title = plot.title,
    column_title_gp = gpar(fontsize=16),
    annotation_legend_list = lgd,
    annotation_legend_side = "left"
  )
  dev.off()
}


# Reformat data
cutie_df = filter(cutie_df, grepl('^WB_|^HC_', cutie_df$var2))
lefse_df = filter(lefse_df, p<0.05)
lefse_df$class = mapvalues(lefse_df$class,
                           c('PsychoticCase', 'HealthyControl'),
                           c('Psychotic Case', 'Healthy Control'),
                           warn_missing=F)

# Create annotation label and color maps
imaging_vars = cutie_df$var2
annotmap_top = ifelse(grepl('^WB_|vol', imaging_vars), 'Volume', 'Metabolite') %>%
  setNames(imaging_vars)
annotcol_top = c(Volume='orange', Metabolite='navy')
annotcol_right = c("Psychotic Case"='#dd443c', 
                   "Healthy Control"='#0073C2FF')


# Create heatmap input matrices
cutie_matrices = create_cutie_matrices(cutie_df)
hm_mat = create_heatmap_matrices(cutie.matrices=cutie_matrices,
                                 lefse.filt=lefse_df$feature)

# Create row/col name maps
M_rownames = sapply(rownames(hm_mat$r), rename_taxa)
M_colnames = gsub('_', ' ', colnames(hm_mat$r))

# Plot and save heatmap
fp = glue('heatmap.png')
plot_heatmap(
  hm.inputs=hm_mat, lefse.df=lefse_df,
  fp=fp, f.w=8.5, f.h=7, plot.title=glue("Plot Title"),
  hm.w=4, hm.h=12, hm.title='Heatmap Title', 
  annot.r.lefse.ymax=6,
  hm.rownames=M_rownames,
  hm.colnames=M_colnames,
  row.km=2, col.km=1,
  pch.size=5
)















