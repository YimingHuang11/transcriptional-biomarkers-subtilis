source("data.R")
source("util.R")

#### data patterns #####
filedir <- 'results/DataPatterns'
if (!dir.exists(filedir)) {dir.create(filedir,recursive=TRUE)}

# boxplot
cat('plotting boxplot\n')
pdf(paste0(filedir,'/boxplot.pdf'),width=12, height=8)
boxplot(m_express,xlab='sample',ylab='gene expression level')
dev.off()

# IQR-median plot
cat('plotting IQR-median plot\n')source("util.R")

pdf(paste0(filedir,'/IQR_median.pdf'),width=12, height=8)
topvariant <- df_genes$IQRexp>quantile(df_genes$IQRexp,0.99)
plot(df_genes$medianexp,df_genes$IQRexp,pch = 19,cex=0.5,main='gene IQR-meadian plot',xlab='median',ylab='IQR')
text(x = df_genes$medianexp[topvariant]-0.1,y = df_genes$IQRexp[topvariant],labels = df_genes$Name[topvariant],cex = 0.5,col='blue')
dev.off()

# gene expression distribution across samples
cat('plotting distribution density plot\n')
plotdistribution<-function(data,gene_list,firedir){
  xmin=min(data)
  xmax=max(data)
  n<-length(gene_list)
  for (i in 1:n){
    pdf(paste0(filedir,'/distribution_',gene_list[i],'.pdf'),width=12, height=8)
    express_density<-density(data[gene_list[i],])
    plot(express_density,xlim=c(xmin,xmax),ylim=c(0,2.2),main=gene_list[i],xlab='Gene expression level', ylab='Density',
         cex.main=3,cex.lab=2,cex.axis=1.5)
    dev.off()

  }
}
plotdistribution(m_express,sample(gene_names,5),firedir)
plotdistribution(m_express,df_genes$Name[order(df_genes$IQRexp)[1:5]],firedir)
plotdistribution(m_express,df_genes$Name[order(df_genes$IQRexp,decreasing=TRUE)[1:5]],firedir)


#########################find differentially expressed genes for condition contrasts##########################################
###conditions contrasts including treatment conditions versus control conditions and similar conditions in different experiments

# matrix of design matrix  - num_samples*num_conditions: '1' treatment '0' control 'non' irrelevant to this condition
m_index<-read.table('config/design_matrix.csv', sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
condition_contrasts<-colnames(m_index)

# run limma modified t-tests for all condition contrasts
# output matrix of log fold changes - num_genes*num_conditions, p-values - num_genes*num_conditions, regulation signs - num_genes*num_conditions: 1 upregualte, -1 downregulate, 0 non-significant
filedirec='results/DEGs'
if (!dir.exists(filedirec)) {dir.create(filedirec,recursive=TRUE)}
DEG_patterns(m_express,m_index,block_list=NULL,condition_list=condition_contrasts,filedirec=filedirec,p.value=0.05,lfc=1)
m_lfcs<-read.table(paste0(filedirec,'/logfoldchanges.csv'), sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
m_pvaules<-read.table(paste0(filedirec,'/pvalues.csv'), sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
m_signs<-read.table(paste0(filedirec,'/expression_signs.csv'), sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ranges of sign p-values
m_neglog_pvaules<--log10(m_pvaules)
m_sign_pvaules<-m_pvaules
m_sign_pvaules[m_sign_pvaules<0.01]<-2
m_sign_pvaules[m_sign_pvaules<0.05]<-1
m_sign_pvaules[m_sign_pvaules<1]<-0
m_sign_pvaules<-m_signs*m_sign_pvaules

# heatmaps of sign p-values and log fold changes
x<-pheatmap::pheatmap(t(m_lfcs), fontsize_col = 0.5,show_colnames = FALSE,
         file = paste0(filedirec,"/log2_fold_changes.pdf"),color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
         main = 'log fold changes in gene expressions under different conditions',width=12,height = 8)
pheatmap::pheatmap(t(m_sign_pvaules), fontsize_col = 0.5,cluster_rows=x$tree_row,cluster_cols=x$tree_col,show_colnames = FALSE,
                   file = paste0(filedirec,'/sign_pvaules.pdf'),color = c("#313695","#74ADD1", "#FFFFFF", "#F46D43","#A50026"),
                   main = '-sign p-values for genes differentially under different conditions ',width=12,height = 8)
