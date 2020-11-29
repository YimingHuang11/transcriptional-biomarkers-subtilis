# set default values for the arguments
threshold_highp<-0.7
threshold_lowp<-0.3
threshold_varp<-0.3
condition_toremove<-'sporulation late stage'

args = commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  condition_toremove<-args[1]
}
if (length(args)>=3) {
  threshold_highp<-as.numeric(args[1])
  threshold_lowp<-as.numeric(args[2])
  threshold_varp<-as.numeric(args[3])
  if (!(threshold_highp<1)*(threshold_lowp<1)*(threshold_varp<1)*(threshold_highp>0)*(threshold_lowp>0)*(threshold_varp>0)*(threshold_lowp<threshold_highp))
    stop("percentile thresholds setting wrong!", call.=FALSE)
  if(length(args)==4)
    condition_toremove<-args[4]
}

######## Data processing #######
source("data.R")
source("util.R")
### cleaning ###

cat('removing uninterested samples and genes...\n')

## genes highly or lowly expressed in all conditions
getInvariantgenes<-function(data,threshold_highp,threshold_lowp,threshold_varp){
  threshold_highexp<-quantile(data,probs = threshold_highp)
  threshold_lowexp<-quantile(data,probs = threshold_lowp)
  threshold_lowIQR<-quantile(df_genes$IQRexp,probs=threshold_varp)
  # high expressed in all samples while IQR is low
  count_highsammples<-rowSums(data>threshold_highexp)
  allhigh_genes<-gene_names[count_highsammples==num_samples & df_genes$IQRexp<threshold_lowIQR]
  # low expressed in all samples while IQR is low
  count_lowsammples<-rowSums(data<threshold_lowexp)
  alllow_genes<-gene_names[count_lowsammples==num_samples & df_genes$IQRexp<threshold_lowIQR]
  invariantGenes<-union(alllow_genes,allhigh_genes)
  return (invariantGenes)

}
invariantGenes<-getInvariantgenes(data=m_express,threshold_highp=0.7,threshold_lowp=0.3,threshold_varp=0.3)

## genes differentially expressed in the condition to be removed
m_lfcs<-read.table('results/DEGs/logfoldchanges.csv', sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
m_pvaules<-read.table('results/DEGs/pvalues.csv', sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)


if(! condition_toremove %in% colnames(m_pvaules))
  stop("Condition to be removed is incorrect!", call.=FALSE)
genenames_toremove<-gene_names[m_pvaules[,condition_toremove]<0.05 & abs(m_lfcs[,condition_toremove])>1]

## samples under the condition to be removed
if(!condition_toremove %in% sample_list$condition)
  stop("Condition to be removed is incorrect!", call.=FALSE)
SampleID_toremove<-sample_list$SampleID[sample_list$condition==condition_toremove]

## discard invariant genes and condition-to-be-removed-related genes
genenames_reduced<-gene_names[!gene_names%in%union(invariantGenes,genenames_toremove)]
write.table(df_genes[rownames(df_genes)%in%genenames_reduced,c('Name','Locus_tag')],"data/selectedGenes.csv", append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

## discard condition-to-be-removed-related samples
SampleID_reduced<-sample_list$SampleID[!sample_list$SampleID %in% SampleID_toremove]


## reduced gene expression data
m_express_reduced<-m_express[genenames_reduced,SampleID_reduced]
cat("The reduced data set contains", length(genenames_reduced),"genes and",length(SampleID_reduced), "samples\n")

### standardisation  ###

cat('Normalising the data via standardisation with respect to corresponding reference expressions...\n')

# read treatment conditions ID and their corresponding reference samples ID

cat('Reading the treatment samples and reference samples...\n')

df_treatref<-read.table('config/treatment_reference_ID.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
index_treatment<-sample_list$sampleIndex[sample_list$ConditionID %in% df_treatref$treatment_condition]
treatment_conditions<-unique(df_treatref$treatment_condition)

# equation: (1) m (treatment counts) = 2^ m(treatment expression levels); m (reference counts) =2^m (average reference expression levels)
#           (2) m (relative counts = m (treatment counts)  – m (reference counts)
#           (3) m (relative expression levels) = sign[m (relative counts)*log2{abs[m (relative counts)]+1}
m_count<-2^m_express
m_count_relative<-matrix()
for (treatment_condition in treatment_conditions){
  treatment_samples<-sample_list$SampleID[sample_list$ConditionID==treatment_condition]
  m_count_treatment<-m_count[,treatment_samples]
  reference_samples<-df_treatref$reference_samples[df_treatref$treatment_condition==treatment_condition]
  m_count_reference<-matrix(rep(rowMeans(m_count[,reference_samples]),length(treatment_samples)),ncol=length(treatment_samples),
                            byrow = FALSE,dimnames=list(gene_names,treatment_samples))
  if(is.na(m_count_relative))
    m_count_relative<-m_count_treatment-m_count_reference
  else
    m_count_relative<-cbind(m_count_relative,m_count_treatment-m_count_reference)
}
m_count_relative<-m_count_relative[,sample_list$SampleID[sample_list$SampleID%in% colnames(m_count_relative)]]
m_express_rel<-sign(m_count_relative)*log2(abs(m_count_relative)+1)
write.table(m_express_rel, 'data/m_express_rel.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

## processed gene expression data with data cleaning and standardisation
m_express_processed<-m_express_rel[genenames_reduced,]

cat("The processed data set contains", length(genenames_reduced),"genes and",length(treatment_samples), "samples\n")

######### comparison between UMAP embedding of transcriptomics data before and after data processing ############
filedir<-'results/umapCompare'
if (!dir.exists(filedir)) {dir.create(filedir)}
setwd(filedir)

cat('Saving the UMAP embeddings transformed from original data and processed data...\n')

## original data
umap_embedding<- uwot::umap(t(m_express),n_components = 2,min_dist = 0.1,n_neighbors = 15,fast_sgd = FALSE,n_threads=1,verbose=T)
df_umap_original<-as.data.frame(umap_embedding)
colnames(df_umap_original)<-c('UMAP1','UMAP2')
df_umap_original<-cbind(df_umap_original,sample_list)

p<-ggplot(df_umap_original,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=condition,shape=type))+
  theme(legend.title =element_text(size=12),legend.text =element_text(size=8),axis.text =element_text(size=16),axis.title =element_text(size=16))
ggsave(file="umap_original.pdf", plot = p,units="in", width=15, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_original.html")

p<-ggplot(df_umap_original,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=medium,shape=type))+
  theme(legend.title =element_text(size=30),legend.text =element_text(size=22),axis.text =element_text(size=16),axis.title =element_text(size=16))
ggsave(file="umap_original_medium.pdf", plot = p,units="in", width=14, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_original_medium.html")

p<-ggplot(df_umap_original,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=experiment,shape=type))+
  theme(legend.title =element_text(size=8),legend.text =element_text(size=8),axis.text =element_text(size=8),axis.title =element_text(size=8))
ggsave(file="umap_original_experiment.pdf", plot = p,units="in", width=17, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_original_experiment.html")

## reduced data
umap_embedding<- uwot::umap(t(m_express_reduced),n_components = 2,min_dist = 0.1,n_neighbors = 15,pca = 100,pca_center = TRUE,scale=TRUE,
                            fast_sgd = FALSE,n_threads=1,verbose=T)
write.table(umap_embedding, 'umap_embedding_reduced.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
df_umap_reduced<-as.data.frame(umap_embedding)
colnames(df_umap_reduced)<-c('UMAP1','UMAP2')
df_umap_reduced<-cbind(df_umap_reduced,sample_list[sample_list$SampleID %in% SampleID_reduced,])

p<-ggplot(df_umap_reduced,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=condition,shape=type))+
  theme(legend.title =element_text(size=12),legend.text =element_text(size=8),axis.text =element_text(size=16),axis.title =element_text(size=16))
ggsave(file="umap_reduced.pdf", plot = p,units="in", width=15, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_reduced.html")

## processed Data
umap_embedding<- uwot::umap(t(m_express_processed),n_components = 2,min_dist = 0.1,n_neighbors = 15,pca = 100,pca_center = TRUE,scale=TRUE,
                            fast_sgd = FALSE,n_threads=1,verbose=T)
write.table(umap_embedding, 'umap_embedding_processed.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)

df_umap_processed<-as.data.frame(umap_embedding)
colnames(df_umap_processed)<-c('UMAP1','UMAP2')
df_umap_processed<-cbind(df_umap_processed,sample_list[index_treatment,])

p<-ggplot(df_umap_processed,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=condition),shape=17)+
  theme(legend.title =element_text(size=12),legend.text =element_text(size=8),axis.text =element_text(size=16),axis.title =element_text(size=16))
ggsave(file="umap_processed.pdf", plot = p,units="in", width=12, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_processed.html")

p<-ggplot(df_umap_processed,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=medium,shape=type))+
  theme(legend.title =element_text(size=30),legend.text =element_text(size=22),axis.text =element_text(size=16),axis.title =element_text(size=16))
ggsave(file="umap_processed_medium.pdf", plot = p,units="in", width=14, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_processed_medium.html")

p<-ggplot(df_umap_processed,aes(x=UMAP1,y=UMAP2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+
  geom_point(size=3,aes(colour=experiment,shape=type))+
  theme(legend.title =element_text(size=8),legend.text =element_text(size=8),axis.text =element_text(size=8),axis.title =element_text(size=8))
ggsave(file="umap_processed_experiment.pdf", plot = p,units="in", width=17, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_processed_experiment.html")

setwd('../../')

######### construct the transcriptional landscape with processed data ############
filedir<-'results/landscape'
if (!dir.exists(filedir)) {dir.create(filedir)}

pd <- new("AnnotatedDataFrame",data=sample_list[index_treatment,])
fd <- new("AnnotatedDataFrame",data=df_genes[genenames_reduced,])

## dimension reduction with umap
# set 'fast_sgd = FALSE'  'cores = 1' to get the same output for each run
umap_embedding<- uwot::umap(t(m_express_processed),n_components = 2,min_dist = 0.1,n_neighbors = 15,
                             fast_sgd = FALSE,n_threads=1,verbose=T)
write.table(umap_embedding, 'results/landscape/umap_embedding.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = T)
df_umap<-umap_visual(umap_embedding,sample_list=sample_list[index_treatment,],filedir='results/landscape',filetoken="umap_condition")

## leiden clustering
## scan on parameters and cluster the samples in the landscape
cluster_parameters<-read.table('config/cluster_parameters.csv',sep=',',header=TRUE,stringsAsFactors = FALSE)
for (i in 1:dim(cluster_parameters)[1]){
  cluster_lei <- monocle3:::leiden_clustering(umap_embedding, pd,k = cluster_parameters[i,'k'],resolution_parameter = cluster_parameters[i,'resolution_parameter'],num_iter = 3,random_seed = 0,weight=F,verbose = T)
  cluster_identity<-cluster_lei$optim_res$membership
  names(cluster_identity)<-sample_list$SampleID[index_treatment]
  write.table(cluster_identity, paste0('results/landscape/cluster_identity_solution',i,'.csv'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = T)
  df_umap<-umap_visual(umap_embedding,cluster=cluster_identity,sample_list=sample_list[index_treatment,],filedir='results/landscape',filetoken=paste0("umap_cluster_solution",i))
}
