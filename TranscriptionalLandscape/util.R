###### global variables #######
# colours
custom.col <-c('#fabed4','#f58231','#FF2F09', '#dcbeff','#f032e6','#911eb4','#ffe119', '#808000','#000075','#0CF091','#a9a9a9', '#9a6324', '#800000','#33614E', '#4363d8','#bcf60c','#aaffc3','#42d4f4','#68A864','#000000')

## gene annotations exported from subtiwiki combined with
## Antisense information (TableS11 in http://genome.jouy.inra.fr/basysbio/bsubtranscriptome/)
## and regulatory network information (https://www.frontiersin.org/articles/10.3389/fmicb.2016.00275/full)
gene_annotations<-read.csv("data/GeneAnnotations.csv", sep="\t",header = TRUE,quote='',stringsAsFactors = FALSE)
rownames(gene_annotations)<-gene_annotations$gene_name

## regulon list
df_regulons<-read.csv("data/regulons.csv", sep=",",header = TRUE,row.names=NULL,quote='',stringsAsFactors = FALSE)



###### functions ########
# run DE test for the treatment group against control group, cutoff at 'p.value' and 'lfc'
# 'treatment' - treatment group index, 'blocklist' - biological replicate list (NULL then do not remove correlation between duplicate spots. )
# print volcano plot (show names of top 'highlight' p-value genes), MDplot plot (hightlight significant genes if 'highlight' !=0)
# save the upregulated gene list and downregulated gene list, return the top table of significant genes []
DEgenes<-function(data,treatment,block_list=NULL,filedirec,p.value=0.05,lfc=0,highlight=50){

  design <- model.matrix(~treatment)
  fit <- lmFit(data, design)
  if(length(block_list)==0){
    fit <- lmFit(data,design)

  }
  else{
    dupcor <- duplicateCorrelation(data,design,block=block_list)
    fit <- lmFit(data,design,block=block_list,correlation=dupcor$consensus.correlation)
  }
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)
  n<-dim(data)[1]
  diff<-topTable(fit, coef="treatment1",n=n, p.value=p.value,lfc=lfc)

  uplist<-rownames(diff[diff$logFC>0,])
  downlist<-rownames(diff[diff$logFC<0,])
  write.table(uplist, paste0(filedirec,'_uplist.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
  write.table(downlist, paste0(filedirec,'_downlist.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
  write.table(diff, paste0(filedirec,'_DEGs.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

  pdf(paste0(filedirec,'_volcano.pdf'))
  x<-unlist(strsplit(filedirec,split ='/'))
  titlename<-x[length(x)]
  if(highlight)
    volcanoplot(fit,coef = 'treatment1',highlight = highlight,names=rownames(fit$t),main=titlename)
  else
    volcanoplot(fit,coef = 'treatment1',main=titlename)
  abline(v=1, col="blue")
  abline(v=-1, col="blue")
  abline(h=-log10(0.05), col="blue")
  dev.off()

  results <- decideTests(fit,p.value = p.value, lfc=lfc)
  pdf(paste0(filedirec,'_MDplot.pdf'), width = 10 ,height = 8)
  plotMD(fit,coef= 'treatment1',status=results[,2],values=c(1,-1),hl.col=c('red','blue'),hl.cex=0.3)
  x<-fit$Amean
  y<-as.matrix(fit$coefficients)[,'treatment1']
  if(highlight){
    text(x = x[uplist],y = y[uplist],labels = uplist,cex = 0.3,col='red')
    text(x = x[downlist],y = y[downlist],labels = downlist,cex = 0.3,col='blue')
  }

  dev.off()
  return (diff)
}

# run DE test for the treatment group against control group, cutoff at 'p.value' and 'lfc'
# 'treatment' - treatment group index, 'blocklist' - biological replicate list (block_list='NULL', correlation=FALSE then do not remove correlation between duplicate spots. )
# print volcanoplot (show names of top 'highlight' p-value genes), MDplot  (hightlight significant genes if 'highlight' !=0)
# save the upregulated gene list and downregulated gene list only look at p-values, return the whole table of DE test
DEtest<-function(data,treatment,block_list=NULL,correlation=FALSE,filedirec,p.value=1,lfc=0,highlight=10){

  design <- model.matrix(~treatment)
  fit <- lmFit(data, design)
  if(correlation){
    dupcor <- duplicateCorrelation(data,design,block=block_list)
    fit <- lmFit(data,design,block=block_list,correlation=dupcor$consensus.correlation)
    cat('remove correlation between duplicates\n')
  }
  else
    fit <- lmFit(data,design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)
  n<-dim(data)[1]
  diff<-topTable(fit, coef="treatment1",n=n, p.value=p.value,lfc=lfc)

  uplist<-rownames(diff[diff$logFC>0,])
  downlist<-rownames(diff[diff$logFC<0,])
  write.table(uplist, paste0(filedirec,'_uplist.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
  write.table(downlist, paste0(filedirec,'_downlist.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
  if(length(uplist)>0){
    uplist_anno<- cbind(gene_annotations[uplist,c(2:6,8:13)],diff[uplist,c(1,5)])
    write.table(uplist_anno, paste0(filedirec,'_uplist.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)
  }
  if(length(downlist)>0){
    downlist_anno<- cbind(gene_annotations[downlist,c(2:6,8:13)],diff[downlist,c(1,5)])
    write.table(downlist_anno, paste0(filedirec,'_downlist.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)
  }
  write.table(diff, paste0(filedirec,'_DEGs.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

  results <- decideTests(fit,p.value = p.value, lfc=lfc)
  pdf(paste0(filedirec,'_MDplot.pdf'), width = 10 ,height = 8)
  plotMD(fit,coef= 'treatment1',status=results[,2],values=c(1,-1),hl.col=c('red','blue'),hl.cex=0.3)
  x<-fit$Amean
  y<-as.matrix(fit$coefficients)[,'treatment1']
  if(highlight){
    if(length(uplist)>0)
      text(x = x[uplist],y = y[uplist],labels = uplist,cex = 0.5,col='red')
    if(length(downlist)>0)
      text(x = x[downlist],y = y[downlist],labels = downlist,cex = 0.5,col='blue')
  }
  dev.off()

  pdf(paste0(filedirec,'_volcano.pdf'))
  wholetable<-topTable(fit, coef="treatment1",n=n, p.value=1,lfc=0)
  x<-unlist(strsplit(filedirec,split ='/'))
  titlename<-x[length(x)]
  x<-wholetable$logFC
  y<--log10(wholetable$adj.P.Val)
  plot(x,y,pch=16,cex=0.5,xlab = "Log2 Fold Change",ylab='-log10(p-value)',main=titlename,col='grey')
  diff<-topTable(fit, coef="treatment1",n=n, p.value=p.value,lfc=lfc) #  look at p-value and log fold change
  if(length(diff)){
  x<-diff$logFC
  y<--log10(diff$adj.P.Val)
  points(x,y,pch=16,cex=0.5)
  if(highlight){
    uplist_x<-x[x>0]
    uplist_y<-y[x>0]
    uplist_label<-rownames(diff)[x>0]
    downlist_x<-x[x<0]
    downlist_y<-y[x<0]
    downlist_label<-rownames(diff)[x<0]
    if(length(uplist_x)>0)
      text(x = uplist_x[1:highlight],y = uplist_y[1:highlight],labels = uplist_label[1:highlight],cex = 1,col='red')
    if(length(downlist_x)>0)
      text(x = downlist_x[1:highlight],y = downlist_y[1:highlight],labels = downlist_label[1:highlight],cex = 1,col='blue')
  }
  dev.off()

  }
  return (wholetable)
}

# m_index - design matrix(num_samples, num_conditions)
# write matrix (num_samples, num_conditions) of expression upregulated/downregulated signs, fold changes, p-values.
DEG_patterns<-function(data,m_index,block_list=NULL,condition_list,filedirec,p.value=1,lfc=0){

  m_sign<-matrix(0, nrow=length(rownames(data)),ncol=length(condition_list),
                 dimnames = list(rownames(data),condition_list))
  m_fc<-matrix(0, nrow=length(rownames(data)),ncol=length(condition_list),
               dimnames = list(rownames(data),condition_list))
  m_pvalue<-matrix(0, nrow=length(rownames(data)),ncol=length(condition_list),
               dimnames = list(rownames(data),condition_list))
  cat('Differential expression analysis performing for conditions:\n')
  for (i in 1:length(condition_list)){
    condition<-condition_list[i]
    cat(i,condition,'\n')
    index<-m_index[,condition]
    index<-which(index!='non')
    treatment<-as.factor(m_index[index,condition])
    wholetable<-DEtest(data[,index],treatment,filedire=paste0(filedirec,'/',condition),p.value=0.05,lfc=1)
    m_sign[rownames(wholetable[wholetable$logFC>lfc  & wholetable$adj.P.Val<p.value,]),condition]<-1
    m_sign[rownames(wholetable[wholetable$logFC<(-lfc) & wholetable$adj.P.Val<p.value,]),condition]<--1
    m_fc[rownames(wholetable),condition]<-wholetable$logFC
    m_pvalue[rownames(wholetable),condition]<-wholetable$adj.P.Val
  }
  write.table(m_sign,paste0(filedirec,'/expression_signs.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_fc,paste0(filedirec,'/logfoldchanges.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_pvalue,paste0(filedirec,'/pvalues.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
}

DEtest_cluster <- function(data,identity,block_list,clusterID,p.cutoff=0.05,lfc.cutoff=1) {
  identity[identity !=clusterID] <- 0
  identity[identity ==clusterID] <- 1
  ct<-factor(identity)
  design <- model.matrix(~ct)
#  dupcor <- duplicateCorrelation(data,design,block=block_list)
#  fit <- lmFit(data,design,block=block_list,correlation=dupcor$consensus.correlation)
  fit <- lmFit(data,design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)

  n<-dim(data)[1]
  diff<-topTable(fit, coef="ct1",n=n, p.value=p.cutoff,lfc=lfc.cutoff)
  uplist<-rownames(diff[diff$logFC>0,])
  downlist<-rownames(diff[diff$logFC<0,])

  results <- decideTests(fit,p.value = p.cutoff, lfc=lfc.cutoff)
  plotMD(fit,coef= 'ct1',status=results[,2],values=c(1,-1),hl.col=c('red','blue'),hl.cex=0.2)
  x<-fit$Amean
  y<-as.matrix(fit$coefficients)[,'ct1']
  if(length(uplist)>0)
    text(x = x[uplist],y = y[uplist],labels = uplist,cex = 0.5,col='red')
  if(length(downlist)>0)
    text(x = x[downlist],y = y[downlist],labels = downlist,cex = 0.5,col='blue')
  wholetable<-topTable(fit, coef="ct1",n=n, p.value=1,lfc=0)
  return (wholetable)
}

DEG_clusters<-function(data,identity,block_list,p.cutoff=0.05,lfc.cutoff=0,filedirec){

  num_clusters<-length(unique(identity))

  m_sign<-matrix(0, nrow=length(rownames(data)),ncol=num_clusters,
                 dimnames = list(rownames(data),1:num_clusters))
  m_fc<-matrix(0, nrow=length(rownames(data)),ncol=num_clusters,
               dimnames = list(rownames(data),1:num_clusters))
  m_pvalue<-matrix(0, nrow=length(rownames(data)),ncol=num_clusters,
                   dimnames = list(rownames(data),1:num_clusters))
  cat('Differential expression analysis performing for each cluster VS others:\n')

  DE_list<-list()
  for(i in 1:num_clusters){
    pdf(paste0(filedirec,'/MD_cluster',i,'.pdf'), width = 10 ,height = 8)
    wholetable<-DEtest_cluster(data,identity,block_list,clusterID=i,p.cutoff,lfc.cutoff)
    dev.off()
    cat('cluster ',i,'\n')
    uplist<-rownames(wholetable[wholetable$logFC>lfc.cutoff  & wholetable$adj.P.Val<p.cutoff,])
    downlist<-rownames(wholetable[-wholetable$logFC>lfc.cutoff  & wholetable$adj.P.Val<p.cutoff,])
    write.table(uplist
                , paste0(filedirec,'/UP_cluster',i,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
    write.table(downlist
                , paste0(filedirec,'/DOWN_cluster',i,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)

    m_sign[uplist,i]<-1
    m_sign[downlist,i]<--1
    m_fc[rownames(wholetable),i]<-wholetable$logFC
    m_pvalue[rownames(wholetable),i]<-wholetable$adj.P.Val

    ID<-paste0('cluster',i)
    DE_list[[paste0('up_',ID)]]<-uplist
    DE_list[[paste0('down_',ID)]]<-downlist
  }
  write.table(m_sign,paste0(filedirec,'/expression_signs.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_fc,paste0(filedirec,'/fold_changes.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_pvalue,paste0(filedirec,'/pvalues.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

  return(DE_list)
}

DEG_clusterspair <- function(data,identity,block_list,clusterID1,clusterID2,p.cutoff=0.01,lfc.cutoff=1) {
  identity[identity ==clusterID1] <- 0
  identity[identity ==clusterID2] <- 1
  ct<-factor(identity)
  design <- model.matrix(~ct)
  dupcor <- duplicateCorrelation(data,design,block=block_list)
  dupcor$consensus.correlation
  fit <- lmFit(data,design,block=block_list,correlation=dupcor$consensus.correlation)
  #fit <- lmFit(object@data,design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)
  results <- decideTests(fit,p.value = p.cutoff, lfc=lfc.cutoff)
  plotMD(fit,coef= 'ct1',status=results[,2],values=c(1,-1),hl.col=c('red','blue'),hl.cex=0.2)
  x<-fit$Amean
  y<-as.matrix(fit$coefficients)[,'ct1']
  highlight1<-names(results[,2][results[,2]==1])
  if(length(highlight1)>0)
    text(x = x[highlight1],y = y[highlight1],labels = highlight1,cex = 0.5,col='red')
  highlight2<-names(results[,2][results[,2]==-1])
  if(length(highlight2)>0)
    text(x = x[highlight2],y = y[highlight2],labels = highlight2,cex = 0.5,col='blue')
  top<-topTable(fit, coef="ct1",n=5875, p.value=p.cutoff,lfc=lfc.cutoff)
  return (list('fit'=fit,'topgenes'=top))
}



umap_visual<-function (umap_embedding,cluster=NULL,sample_list,filedir,filetoken){

  df_umap<-as.data.frame(umap_embedding)
  colnames(df_umap)<-c('component1','component2')
  df_umap<-cbind(df_umap,sample_list)

  setwd(filedir)

  if (!is.null(cluster)){
    num_clusters<-length(unique(cluster))
    label_clusters<-paste0('cluster_',1:num_clusters)
    df_umap$cl<-factor(paste0('cluster_',as.character(cluster)),levels=label_clusters)
    p<-ggplot(df_umap,aes(x=component1,y=component2,color=cl,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+geom_point(size=2)+scale_color_manual(values = custom.col[1:num_clusters])+theme(legend.title = element_blank())+ggtitle('Clustering samples in umap space')
    ggsave(file=paste0(filetoken,'.pdf'), plot = p, units="in", width=12, height=10, dpi = 300)
    p<-ggplotly(p)
    htmlwidgets::saveWidget(p, file=paste0(filetoken,'.html'))
  }
  else{
  num_conditions<-length(unique(df_umap$condition))
  p<-ggplot(df_umap,aes(x=component1,y=component2,color=condition,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+geom_point(size=2)+theme(legend.title = element_blank())+ggtitle('Samples colored by major conditions in umap space')
  ggsave(file=paste0(filetoken,'.pdf'), plot = p,units="in", width=12, height=10, dpi = 300)
  p<-ggplotly(p)
  htmlwidgets::saveWidget(p, file=paste0(filetoken,'.html'))
  }
  setwd(paste(rep('../',stringr::str_count(filedir, "/")+1),collapse=''))
  return (df_umap)
}

DEtest_cluster <- function(data,identity,clusterID,p.cutoff=0.05,lfc.cutoff=1) {
  identity[identity !=clusterID] <- 0
  identity[identity ==clusterID] <- 1
  ct<-factor(identity)
  design <- model.matrix(~ct)
  fit <- lmFit(data,design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)

  n<-dim(data)[1]
  diff<-topTable(fit, coef="ct1",n=n, p.value=p.cutoff,lfc=lfc.cutoff)
  uplist<-rownames(diff[diff$logFC>0,])
  downlist<-rownames(diff[diff$logFC<0,])

  results <- decideTests(fit,p.value = p.cutoff, lfc=lfc.cutoff)
  plotMD(fit,coef= 'ct1',status=results[,2],values=c(1,-1),hl.col=c('red','blue'),hl.cex=0.2)
  x<-fit$Amean
  y<-as.matrix(fit$coefficients)[,'ct1']
  if(length(uplist)>0)
    text(x = x[uplist],y = y[uplist],labels = uplist,cex = 0.5,col='red')
  if(length(downlist)>0)
    text(x = x[downlist],y = y[downlist],labels = downlist,cex = 0.5,col='blue')
  wholetable<-topTable(fit, coef="ct1",n=n, p.value=1,lfc=0)
  return (wholetable)
}
DEG_clusters<-function(data,identity,p.cutoff=0.05,lfc.cutoff=0,filedirec){

  dir.create(filedirec, showWarnings = FALSE)
  num_clusters<-length(unique(identity))
  m_sign<-matrix(0, nrow=length(rownames(data)),ncol=num_clusters,dimnames = list(rownames(data),1:num_clusters))
  m_fc<-matrix(0, nrow=length(rownames(data)),ncol=num_clusters,dimnames = list(rownames(data),1:num_clusters))
  m_pvalue<-matrix(0, nrow=length(rownames(data)),ncol=num_clusters,dimnames = list(rownames(data),1:num_clusters))
  DE_list<-list()
  for(i in 1:num_clusters){
    pdf(paste0(filedirec,'/MD_cluster',i,'.pdf'), width = 10 ,height = 8)
    wholetable<-DEtest_cluster(data,identity,clusterID=i,p.cutoff,lfc.cutoff)
    dev.off()
    uplist<-rownames(wholetable[wholetable$logFC>lfc.cutoff  & wholetable$adj.P.Val<p.cutoff,])
    downlist<-rownames(wholetable[-wholetable$logFC>lfc.cutoff  & wholetable$adj.P.Val<p.cutoff,])
    write.table(uplist
                , paste0(filedirec,'/UP_cluster',i,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
    write.table(downlist
                , paste0(filedirec,'/DOWN_cluster',i,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)

    m_sign[uplist,i]<-1
    m_sign[downlist,i]<--1
    m_fc[rownames(wholetable),i]<-wholetable$logFC
    m_pvalue[rownames(wholetable),i]<-wholetable$adj.P.Val

    ID<-paste0('cluster',i)
    DE_list[[paste0('up_',ID)]]<-uplist
    DE_list[[paste0('down_',ID)]]<-downlist
  }
  write.table(m_sign,paste0(filedirec,'/expression_signs.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_fc,paste0(filedirec,'/fold_changes.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_pvalue,paste0(filedirec,'/pvalues.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

  return(DE_list)
}
ARG_clusters<-function(data,identity,upthreshold,downthreshold,percent_threshold,filedirec){
  dir.create(filedirec, showWarnings = FALSE)
  clustersizes<-table(identity)
  n_clusters<-length(unique(identity))
  genenames<-rownames(data)
  AR_list<-list()
  m_percent<-matrix(0, nrow=dim(data)[1],ncol=2*n_clusters,dimnames = list(rownames(data),c(rbind(paste0('up',1:n_clusters), paste0('down',1:n_clusters)))))
  m_avg<--matrix(0, nrow=dim(data)[1],ncol=n_clusters,dimnames = list(rownames(data),paste0('cluster',1:n_clusters)))
  for (clusterno in 1:n_clusters){
    current_uppercent<-rowSums(data[,which(identity==clusterno)]>upthreshold)/clustersizes[clusterno]
    uplist<-genenames[current_uppercent>=percent_threshold]
    write.table(uplist,paste0(filedirec,'/UP_cluster',clusterno,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
    m_percent[,clusterno*2-1]<-current_uppercent

    current_downpercent<-rowSums(data[,which(identity==clusterno)]<downthreshold)/clustersizes[clusterno]
    downlist<-genenames[current_downpercent>=percent_threshold]
    write.table(downlist,paste0(filedirec,'/DOWN_cluster',clusterno,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
    m_percent[,clusterno*2]<-current_downpercent

    m_avg[,clusterno]<-rowMeans(data[,which(identity==clusterno)])

    ID<-paste0('cluster',clusterno)
    AR_list[[paste0('up_',ID)]]<-uplist
    AR_list[[paste0('down_',ID)]]<-downlist
  }
  write.table(m_percent,paste0(filedirec,'/percentages.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_avg,paste0(filedirec,'/intensities.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

  return (AR_list)

}

write.markeranno<-function(maker_lists,p_values=NULL,lfcs=NULL,percentages=NULL,intensities=NULL,filename,n=10){

  df_makerannos <- data.frame(clusterID=numeric(),mode=character(),gene_name=character(),
                              p_value=numeric(),fold_change=numeric(),percentage=numeric(),intensity=numeric(),stringsAsFactors=FALSE)
  clusterIDs<-names(maker_lists) #both up and down, length(clusterIDs)=2*num_clusters
  modes<-c('upregulated','downregulated')
  for (i in 1:length(clusterIDs)){
    temp_genes<-maker_lists[[clusterIDs[i]]]
    if(length(temp_genes)>=n)
      temp_genes<-temp_genes[1:n]
    temp_df<-data.frame(rep(ceiling(i/2),length(temp_genes)),rep(modes[(i+1)%%2+1],length(temp_genes)),temp_genes,
                        p_values[temp_genes,ceiling(i/2)],lfcs[temp_genes,ceiling(i/2)],percentages[temp_genes,i],intensities[temp_genes,ceiling(i/2)])
    colnames(temp_df)<-c('clusterID','mode','gene_name','p_value','fold_change','percentage','intensity')
    df_makerannos = rbind(df_makerannos,temp_df)
  }
  df_makerannos<-cbind(df_makerannos,gene_annotations[match(df_makerannos$gene_name,gene_annotations$gene_name),3:13])
  write.table(df_makerannos,filename, append = FALSE, sep = "\t",quote=FALSE,col.names = T, row.names = F)
  return(df_makerannos)
}
write.sampleanno<-function(cluster,sample_list,filename){
  df_annos <- data.frame(clusterID=numeric(),sample_annotation=character(),majoy_treatment=character(),stringsAsFactors=FALSE)
  for (i in 1:length(unique(cluster))){
    indices<-as.numeric(names(cluster)[which(cluster==i)])
    annotations<-unique(sample_list$annotation[indices])
    temp_df<-data.frame(rep(i,length(annotations)),annotations,sample_list$condition[match(annotations,sample_list$annotation)])
    names(temp_df)<-c('clusterID','sample_annotation','majoy_treatment')
    df_annos = rbind(df_annos,temp_df)
  }
  write.table(df_annos,filename, append = FALSE, sep = "\t",quote=FALSE,col.names = T, row.names = F)
}

Regulator_clusters<-function(data,cluster,df_regulators,filetoken){

  m_regexp<-matrix(0, nrow=length(cluster),ncol=length(unique(df_regulators$regulator)),
                 dimnames = list(names(cluster),sort(unique(df_regulators$regulator))))

  m_percentages<-matrix(0, nrow=length(unique(cluster)),ncol=length(unique(df_regulators$regulator)),
                        dimnames = list(1:length(unique(cluster)),sort(unique(df_regulators$regulator))))
  m_intensities<-matrix(0, nrow=length(unique(cluster)),ncol=length(unique(df_regulators$regulator)),
                        dimnames = list(1:length(unique(cluster)),sort(unique(df_regulators$regulator))))

  genenames<-rownames(data)

  top30_chromosome<-apply(data,2,quantile,probs=0.7)
  bottom30_chromosome<-apply(data,2,quantile,probs=0.3)
  for(c in 1:dim(m_intensities)[2]){

    regulated_genes<-df_regulators$gene[df_regulators$regulator==colnames(m_intensities)[c]]
    regulated_genes<-regulated_genes[regulated_genes%in%genenames]

    if (length(regulated_genes)>1)
      m_regexp[,c]<-colMeans(data[regulated_genes,])
    else
      m_regexp[,c]<-data[regulated_genes,]

    for(r in 1:dim(m_intensities)[1]){
      sample_indices<-names(which(cluster==r))
      upthreshold_exp<-top30_chromosome[sample_indices]
      downthreshold_exp<-bottom30_chromosome[sample_indices]
      if(length(regulated_genes)>1){
        # mean expression in cluster r and regulator c
        exp_regcluster<-data[regulated_genes,sample_indices]
        m_intensities[r,c]<-mean(exp_regcluster) # mean expression in cluster r and regulator c
        # percentage of samples whose regulator mean expression higher than top 30 expressions in chromosome
        regulator_exp<-colMeans(exp_regcluster)
        if(m_intensities[r,c]>=0)
          m_percentages[r,c]<-sum(regulator_exp>=upthreshold_exp)/length(sample_indices)
        else
          m_percentages[r,c]<-sum(regulator_exp<=downthreshold_exp)/length(sample_indices)

        }
    }
  }
  write.table(m_intensities,paste0(filetoken,'intensities.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
  write.table(m_percentages,paste0(filetoken,'percentages.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

  return (m_regexp)
}
