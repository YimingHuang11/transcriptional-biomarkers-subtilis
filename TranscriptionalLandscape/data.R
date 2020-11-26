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

source("util.R")

#### read data #####

# read log2 transformed and quantile normalised gene expression matrix (TableS1 in http://genome.jouy.inra.fr/basysbio/bsubtranscriptome/)
df_geneExpre<-read.table("data/GeneExpre_norm.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
df_geneExpre <- df_geneExpre[,-1]
m_express<-t(as.matrix(df_geneExpre))
num_genes<-dim(m_express)[1]
num_samples<-dim(m_express)[2]
print("Reading gene expression data...")

# read sample annotations
sample_list<-read.table('data/sample_list.csv', sep = ",", header = TRUE,stringsAsFactors = FALSE)
sample_list$experiment<-factor(sample_list$experiment,levels=unique(sample_list$experiment))
colnames(m_express)<-sample_list$SampleID

# read gene annotations
gene_names<-rownames(m_express)
df_genes<-read.table("data/gene_list.csv", sep=",",quote = "", header = TRUE, stringsAsFactors = FALSE)
rownames(df_genes)<-gene_names
df_genes$meanexp<-rowMeans(m_express)
df_genes$medianexp<-matrixStats::rowMedians(m_express)
df_genes$IQRexp<-matrixStats::rowIQRs(m_express)
df_genes$skewness<-apply(m_express,1,e1071::skewness)
df_genes$range<-apply(m_express,1,max)-apply(m_express,1,min)

## gene name and locus tag converting table in Nicolus data
gene_list<-read.table("data/locus_symbol_Nicolus.csv",quote="", sep=",",header = TRUE, stringsAsFactors = FALSE)

## gene annotations exported from subtiwiki combined with
## Antisense information (TableS11 in http://genome.jouy.inra.fr/basysbio/bsubtranscriptome/)
## and regulatory network information (https://www.frontiersin.org/articles/10.3389/fmicb.2016.00275/full)
gene_annotations<-read.csv("data/GeneAnnotations.csv", sep="\t",header = TRUE,quote='',stringsAsFactors = FALSE)
rownames(gene_annotations)<-gene_annotations$gene_name


cat("This data set contains", num_genes,"genes and",num_samples, "samples\n")
cat('Including',length(unique(sample_list$annotation)), "conditions which were conducted in",length(unique(sample_list$experiment)), "experiments by",length(unique(sample_list$laboratory)), "labs\n")
