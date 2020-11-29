########packages######
library(ggplot2)
library(plotly)
library(cowplot)
library(RColorBrewer)
library(plotrix)
library(e1071)
library(dplyr)
library(data.table)
library(igraph)
library(pheatmap)
library(heatmaply)
library(ComplexHeatmap)
library(limma)
library(Cairo)
suppressMessages(library(devtools))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
library(uwot)
library(monocle3)



#### read data #####

# read log2 transformed and quantile normalised gene expression matrix (TableS1 in http://genome.jouy.inra.fr/basysbio/bsubtranscriptome/)
df_geneExpre<-read.table("data/GeneExpre.csv", sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
m_express<-as.matrix(df_geneExpre)
num_genes<-dim(m_express)[1]
num_samples<-dim(m_express)[2]
cat("Reading gene expression data...\n")

# read sample list
sample_list<-read.table('data/sample_list.csv', sep = ",", check.names = FALSE, header = TRUE,stringsAsFactors = FALSE)
sample_list$experiment<-factor(sample_list$experiment,levels=unique(sample_list$experiment))

# read gene list
df_genes<-read.table("data/gene_list.csv", sep=",",quote = "", header = TRUE, stringsAsFactors = FALSE)
gene_names<-rownames(df_genes)
df_genes$meanexp<-rowMeans(m_express)
df_genes$medianexp<-matrixStats::rowMedians(m_express)
df_genes$IQRexp<-matrixStats::rowIQRs(m_express)
df_genes$skewness<-apply(m_express,1,e1071::skewness)
df_genes$range<-apply(m_express,1,max)-apply(m_express,1,min)


cat("This data set contains", num_genes,"genes and",num_samples, "samples\n")
cat('Including',length(unique(sample_list$annotation)), "conditions which were conducted in",length(unique(sample_list$experiment)), "experiments by",length(unique(sample_list$laboratory)), "labs\n")
