# Plotting for specific features

# load libraries
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)

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


# choose colors
col1 <- c("#929aab", "#ce2525")
# col1unpaired <- c("green4", "purple")
col1unpaired <- c("green4", "lightblue")
col2 <- c("#f3a333", "#0074e4", "#8f8787")

# choose line types
line1 <- c("solid", "dashed", "dotted")

# set output directory
dir="~/Desktop/clemente_lab/Projects/oa/outputs/jobs33/"

# for taxa in Limivivens, Anaerostipes, Akkermansia
# probably create a custom df in python


df=read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/inputs/df_R_boxplots.tsv',
              sep='\t', header=TRUE, row.names = 1)

# rewrite order of factors
df$Timepoint <- factor(df$Timepoint, levels = c('Pre', 'Post'))


# rename
df[df=="No response"]<-"NR"
df[df=="Response"]<-"R"

col1 = c('darkgray','red2') 

col2 <- c("#f3a333", "#0074e4", "#8f8787")

# choose line types
line1 <- c("solid", "dashed", "dotted")


taxas = c(
  'Lachnospiraceae_Limivivens',
  'Lachnospiraceae_Anaerostipes',
  'Lachnospiraceae_Butyribacter',
  'Akkermansiaceae_Akkermansia'
)

panels = list()

for (i in 1:length(taxas)){
  taxa = taxas[i]
  
  # paired boxplots for differences
  # https://stackoverflow.com/questions/76569716/r-how-to-make-a-boxplot-with-lines-connecting-paired-points-between-them
  figure <- ggplot(data = df, aes(fill = Timepoint, x=interaction(Timepoint, WOMAC_P_Response, sep = "_"),
                                  y = .data[[taxa]]))+# y = Lachnospiraceae_Butyribacter))+
    geom_point()+ 
    geom_line(aes(group = HostSubjectId), color = "black",linetype = 2) + 
    #scale_linetype_manual(values = line1, name = "Difference", labels = c("Increase", "Decrease", "No change")) +
    #scale_color_manual(values = col2, name = "Difference", labels = c("Increase", "Decrease", "No change")) +      
    geom_boxplot(alpha=0.2) +
    scale_x_discrete(labels = c('Pre_NR', 'Post_NR', 'Pre_R','Post_R')) +
    scale_fill_manual(NULL, values = col1) +      
    theme_bw(base_size = 18) +
    theme(axis.title.x = element_blank(),
          legend.position = "right",
          legend.title=element_blank())
  
  
  fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs33/',taxa,'.pdf')
  pdf(file = fbp, height = 4, width = 6)
  plot(figure)
  dev.off()

  panels[[i]] <- figure
}


# scatterplots
pairs = list(
  c('Gut_Microbiome_Alpha_Diversity','WOMAC_Pain'),
  c('Plasma_Metabolite_Alpha_Diversity','WOMAC_Pain'),
  c('Gut_Lachnospiraceae_Limivivens','WOMAC_Pain')
)


df=read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/inputs/df_corr_R_label.tsv',
              sep='\t', header=TRUE, row.names = 1)
df[df=="No response"]<-"NR"
df[df=="Response"]<-"R"

# rewrite order of factors
# df$Timepoint <- factor(df$Timepoint, levels = c('Pre', 'Post'))

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

for (i in 1:length(pairs)){
  p1 = pairs[[i]][1]
  p2 = pairs[[i]][2]
  fit = lm(df[[p2]] ~ df[[p1]])
  r2=round(summary(fit)$r.squared,2)
  p=round(summary(fit)$coefficients[,4],4)
  
  figure <- ggplot(data = df, aes(x = .data[[p1]], y = .data[[p2]])) + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
    geom_point(size=4,aes(color = WOMAC_P_Response)) + #shape = sib_02),size=4) +
    scale_color_manual(values = c("NR" = col1unpaired[1], "R" = col1unpaired[2])) + #col1, labels = c("Unaffected", "RA")) +
    guides(shape = "none") + 
    geom_smooth(method='lm',color='black') +#,aes(x=.data[[p2]],y= .data[[p1]])) + 
    # ggtitle(paste0('p=',pv)) +
    scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
    scale_y_continuous(labels = f.dec) +
    bbkg + 
    ggtitle(paste0('r-sq=',as.character(r2),', p=',as.character(p)))   # https://stackoverflow.com/questions/39335005/add-p-value-and-r-on-ggplot-follow-up
    # 2 decimal places on y-axis
  
  # save plot
  fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs33/',p1,'_',p2,'.pdf')
  pdf(file = fbp, height = 4, width = 6)
  plot(figure)
  dev.off()
  
  panels[[i + length(taxas)]] <- figure
}


figure <- plot_grid(panels[[1]],
                    panels[[2]],
                    panels[[3]],
                    panels[[4]],
                    labels = c('A', 'B', 'C', 'D'),
                    ncol=2,nrow=2,
                    rel_heights = c(1,1,1,1),
                    rel_widths = c(1,1,1,1))
fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs33/barplots.pdf')
pdf(file = fbp, height = 12, width = 16)
plot(figure)
dev.off()

figure <- plot_grid(panels[[5]],
                    panels[[6]],
                    panels[[7]],
                    labels = c('E','F','G'),
                    ncol=3,nrow=1,
                    rel_heights = c(1,1,1),
                    rel_widths = c(1,1,1))
fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs33/lmplots.pdf')
pdf(file = fbp, height = 4, width = 16)
plot(figure)
dev.off()

# heatmap

df=read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/inputs/df_corr_R.tsv',
              sep='\t', header=TRUE, row.names = 1)

cormat <- round(cor(df),2)

library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri <- cormat #get_upper_tri(cormat)

library(reshape2)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
#cormat <- reorder_cormat(cormat)
#upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

figure <- ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = 'right',
    #legend.position = c(0.6, 0.7),
    legend.direction = "vertical")+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                               title.position = "top", title.hjust = 0.5))


fbp = paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/outputs/jobs33/heatmap.pdf')
pdf(file = fbp, height = 8, width = 8)
plot(figure)
dev.off()






# GLM funsies
df=read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/oa/inputs/df_curated_R.tsv',
              sep='\t', header=TRUE, row.names = 1)

# factor conversions
df$Sex = as.factor(df$Sex)
df$HTN = as.factor(df$HTN)
df$DM = as.factor(df$DM)
df$Race = as.factor(df$Race)

# DM is all 0's
model=glm(WOMAC_Pain ~ Sex + HTN + Race,data=df)
summary(model)

# Anti_Inflammatory_Diet_Score, HTN, Race only 1 AA
# men tended to respond less well
model=glm(WOMAC_Pain ~  Plasma_Metabolite_Alpha_Diversity + Gut_Microbiome_Alpha_Diversity + Sex,data=df)
summary(model)

model=glm(WOMAC_Pain ~  Plasma_Metabolite_Alpha_Diversity + Gut_Microbiome_Alpha_Diversity,data=df)
summary(model)






