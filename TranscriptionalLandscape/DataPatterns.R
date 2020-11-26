source("data.R")
#### data patterns #####
filedir <- 'results/DataPatterns'
if (!dir.exists(filedir)) {dir.create(filedir,recursive=TRUE)}

# boxplot
print('plotting boxplot')
pdf(paste0(filedir,'/boxplot.pdf'),width=12, height=8)
boxplot(m_express,xlab='sample',ylab='gene expression level')
dev.off()

# IQR-median plot
print('plotting IQR-median plot')
pdf(paste0(filedir,'/IQR_median.pdf'),width=12, height=8)
topvariant <- df_genes$IQRexp>quantile(df_genes$IQRexp,0.99)
plot(df_genes$medianexp,df_genes$IQRexp,pch = 19,cex=0.5,main='gene IQR-meadian plot',xlab='median',ylab='IQR')
text(x = df_genes$medianexp[topvariant]-0.1,y = df_genes$IQRexp[topvariant],labels = df_genes$Name[topvariant],cex = 0.5,col='blue')
dev.off()

# gene expression distribution across samples
print('plotting distribution density plot')
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
