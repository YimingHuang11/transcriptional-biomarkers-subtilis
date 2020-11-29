
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Cluster identity input file name must be supplied", call.=FALSE)
}
clusterfile<-args[1]

source("data.R")
source("util.R")

#### read processed data and landscape cluster solution ####
cat('Reading processed data\n')
m_express_rel<-as.matrix(read.table('data/m_express_rel.csv',sep=",", header = TRUE, check.names = FALSE,stringsAsFactors = FALSE))
genenames_reduced<-rownames(read.table('data/selectedGenes.csv',sep=",", header = TRUE, stringsAsFactors = FALSE))

index_treatment<-sample_list$sampleIndex[sample_list$SampleID %in% colnames(m_express_rel)]
m_express_processed<-m_express_rel[genenames_reduced,]

cat('Reading cluster solution\n')
cluster_identity<-read.table(clusterfile,sep=",", header = FALSE, stringsAsFactors = FALSE)
cluster<-cluster_identity[,2]
names(cluster)<-cluster_identity[,1]
n_clusters<-length(unique(cluster))

filedir<-'results/landscape'

#### identify cluster markers ####
cat('Identifying the cluster markers...\n')

filedirec<-paste0(filedir,"/DEtests")

## genes upregulated/downregulated in each cluster vs others (p-value<0.05, lfc>1)
DElists<-DEG_clusters(data=m_express_rel,identity=cluster,p.cutoff=0.05,lfc.cutoff=1,filedirec = filedirec)
p_values<-as.matrix(read.table(paste0(filedirec,'/pvalues.csv'),sep=',',stringsAsFactors = FALSE,header = T))
lfcs<-as.matrix(read.table(paste0(filedirec,'/fold_changes.csv'),sep=',',stringsAsFactors = FALSE,header = T))

## genes show positive/negative relative expressions in 90% samples in each cluster
filedirec<-paste0(filedir,"/ARtests")
ARlists<-ARG_clusters(data=m_express_rel,identity=cluster,upthreshold=0,downthreshold=0,percent_threshold=0.90,filedirec = filedirec)
percentages<-as.matrix(read.table(paste0(filedirec,'/percentages.csv'),sep=',',stringsAsFactors = FALSE,header = T))
intensities<-as.matrix(read.table(paste0(filedirec,'/intensities.csv'),sep=',',stringsAsFactors = FALSE,header = T))

## cluster markers -intersect of (DElists,ARlists)
filedirec<-paste0(filedir,"/cluster_markers")
dir.create(filedirec, showWarnings = FALSE)
cluster_markers<-list()
for(i in 1:n_clusters){
  cluster_markers[[paste0('up_cluster',i)]]<-intersect(DElists[[paste0('up_cluster',i)]],ARlists[[paste0('up_cluster',i)]])
  cluster_markers[[paste0('down_cluster',i)]]<-intersect(DElists[[paste0('down_cluster',i)]],ARlists[[paste0('down_cluster',i)]])
  write.table(cluster_markers[[paste0('up_cluster',i)]],paste0(filedirec,'/UP_cluster',i,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
  write.table( cluster_markers[[paste0('down_cluster',i)]],paste0(filedirec,'/DOWN_cluster',i,'.txt'), append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
}
marker_full<-write.markeranno(cluster_markers,p_values,lfcs,percentages,intensities,filename= paste0(filedirec,"/markers_annotation.csv"),n=1000)
marker_top<-write.markeranno(cluster_markers,p_values,lfcs,percentages,intensities,filename= paste0(filedirec,"/markersTop200_annotation.csv"),n=200)


############# activities of sigma factors, regulators for each cluster########
cat("Calculating the activity levels of different sigma factors, regulators for each cluster\n")

## read sigma factor list and transcriptional regulator list
df_regulator_neg<-df_regulons[!startsWith(df_regulons$regulator,'Sig') & df_regulons$flag==-1,]
df_regulator_pos<-df_regulons[!startsWith(df_regulons$regulator,'Sig') & df_regulons$flag==1,]
df_sigma<-df_regulons[startsWith(df_regulons$regulator,'Sig'),]

dir.create(paste0(filedir,"/regulators"), showWarnings = FALSE)

m_sigexp<-Regulator_clusters(data=m_express_rel,cluster=cluster,df_regulators=df_sigma,filetoken=paste0(filedir,"/regulators/sigmafactor_"))
m_regposexp<-Regulator_clusters(data=m_express_rel,cluster=cluster,df_regulators=df_regulator_pos,filetoken=paste0(filedir,"/regulators/regpos_"))
colnames(m_regposexp)<-paste0(colnames(m_regposexp),'_pos')
m_regnegexp<-Regulator_clusters(data=m_express_rel,cluster=cluster,df_regulators=df_regulator_neg,filetoken=paste0(filedir,"/regulators/regneg_"))
colnames(m_regnegexp)<-paste0(colnames(m_regnegexp),'_neg')
m_regexp<-cbind(m_regposexp,m_regnegexp)
m_regexp<-m_regexp[,sort(colnames(m_regexp))]
index_bycluster<-names(sort(cluster))
m_sigTF_exp<-t(cbind(m_sigexp[index_bycluster,],m_regexp[index_bycluster,]))
write.table(m_sigTF_exp,paste0(filedir,"/regulators/sigTF_exp.csv"),append = FALSE, sep = ",",quote=FALSE, col.names = TRUE, row.names = TRUE)

cat('Plot heatmaps of regulator activities in different clusters\n')

anno_row<-data.frame(row.names=index_bycluster,cluster=factor(sort(cluster),levels=1:n_clusters))
data<-m_regexp[index_bycluster,sort(colnames(m_regexp))]
pheatmap::pheatmap(data,  cluster_rows = F,cluster_cols=F,show_colnames = T, show_rownames = T,
         annotation_row  =anno_row,angle_col = '90',fontsize_col = 6, fontsize_row = 6,
         file = paste0(filedir,"/regulators/hp_regulators_exp.pdf"),width=16,height=16)
anno_hov<-matrix(paste0('[',rep(rownames(data),each=dim(data)[2]),', ',rep(colnames(data),dim(data)[1]),']  = ', c(t(data))),nrow =dim(data)[1],byrow = TRUE)
heatmaply(data,method='plotly',Rowv=FALSE,Colv=FALSE,row_side_colors=anno_row,
          column_text_angle=90,fontsize_row = 6, colorbar_len=0.1,fontsize_col = 6,xlab='regulators',ylab='samples',
           file = paste0(filedir,"/regulators/hp_regulators_exp.html"))

cat('Plot heatmaps of sigma factor activities in different clusters\n')

m_sigma_percentage<-as.matrix(read.table(paste0(filedir, "/regulators/sigmafactor_percentages.csv"), sep=",", header = TRUE, stringsAsFactors = FALSE))
m_sigma_intensity<-as.matrix(read.table(paste0(filedir, "/regulators/sigmafactor_intensities.csv"), sep=",", header = TRUE, stringsAsFactors = FALSE))
df_plotsigma<-data.frame(cluster=rownames(m_sigma_percentage)[row(m_sigma_percentage)], sigma_factor=colnames(m_sigma_percentage)[col(m_sigma_percentage)],
           percentage=c(m_sigma_percentage))
df_plotsigma$intensity<-c(m_sigma_intensity)
df_plotsigma$cluster<-factor(df_plotsigma$cluster,levels=as.character(n_clusters:1))



m_regpos_percentage<-as.matrix(read.table(paste0(filedir,'/regulators/regpos_percentages.csv'), sep=",", header = TRUE, stringsAsFactors = FALSE))
m_regpos_intensity<-as.matrix(read.table(paste0(filedir,'/regulators/regpos_intensities.csv'), sep=",", header = TRUE, stringsAsFactors = FALSE) )
m_regneg_percentage<-as.matrix(read.table(paste0(filedir,'/regulators/regneg_percentages.csv'), sep=",", header = TRUE, stringsAsFactors = FALSE))
m_regneg_intensity<-as.matrix(read.table(paste0(filedir,'/regulators/regneg_intensities.csv'), sep=",", header = TRUE, stringsAsFactors = FALSE) )

df_plotreg_pos<-data.frame(cluster=rownames(m_regpos_percentage)[row(m_regpos_percentage)], regulator=colnames(m_regpos_percentage)[col(m_regpos_percentage)],
                         percentage=c(m_regpos_percentage))
df_plotreg_pos$intensity<-c(m_regpos_intensity)
df_plotreg_pos$cluster<-factor(df_plotreg_pos$cluster,levels=as.character(1:n_clusters))
df_plotreg_neg<-data.frame(cluster=rownames(m_regneg_percentage)[row(m_regneg_percentage)], regulator=colnames(m_regneg_percentage)[col(m_regneg_percentage)],
                       percentage=c(m_regneg_percentage))
df_plotreg_neg$intensity<-c(m_regneg_intensity)
df_plotreg_neg$cluster<-factor(df_plotreg_neg$cluster,levels=as.character(n_clusters:1))

pdf(paste0(filedir,'/regulators/dotplot_reg_pos.pdf'),width=30, height=6)
ggplot(df_plotreg_pos,mapping = aes(x = regulator, y = cluster) )+ geom_point(mapping = aes(size = percentage, colour = intensity)) +
  scale_colour_gradient2(low = "blue",mid = "white",high = "red")+scale_size(breaks=c(0.2,0.4,0.6,0.8,1),range=c(0,10))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  labs(x = 'Transcriptional Regulator',y = 'cluster')
dev.off()

pdf(paste0(filedir,'/regulators/dotplot_reg_neg.pdf'),width=36, height=6)
ggplot(df_plotreg_neg,mapping = aes(x = regulator, y = cluster) )+ geom_point(mapping = aes(size = percentage, colour = intensity)) +
  scale_colour_gradient2(low = "blue",mid = "white",high = "red")+scale_size(breaks=c(0.2,0.4,0.6,0.8,1),range=c(0,10))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  labs(x = 'Transcriptional Regulator',y = 'cluster')
dev.off()

df_plotreg_neg$regulator<-paste0(df_plotreg_neg$regulator,'_neg')
df_plotreg_pos$regulator<-paste0(df_plotreg_pos$regulator,'_pos')
df_plotreg<-rbind(df_plotreg_pos,df_plotreg_neg)
df_plotreg$cluster<-factor(df_plotreg$cluster,levels=as.character(1:n_clusters))
pdf(paste0(filedir,'/regulators/dotplot_reg.pdf'),width=60, height=6)
ggplot(df_plotreg,mapping = aes(x = regulator, y = cluster))+ geom_point(mapping = aes(size = percentage, colour = intensity)) +
  scale_colour_gradient2(low = "blue",mid = "white",high = "red")+scale_size(breaks=c(0.2,0.4,0.6,0.8,1),range=c(0,10))+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  labs(x = 'Transcriptional Regulator',y = 'cluster')
dev.off()


sigma_list<-unique(df_plotsigma$sigma_factor)
m_pi<-matrix(rep(0,length(sigma_list)*n_clusters*2),ncol=n_clusters*2)
for (i in 1:length(sigma_list)){
  x<-c(df_plotsigma[df_plotsigma$sigma_factor ==sigma_list[i],3:4])
  m_pi[i,]<-c(rbind(x[[1]],x[[2]]))
}
colnames(m_pi)<-paste0('cluster',rep(1:n_clusters,each=2),'_',rep(c('percentage','intensity'),n_clusters))
rownames(m_pi)<-sigma_list
write.table(m_pi,paste0(filedir,'/regulators/sigma_sumamry.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = TRUE, row.names = TRUE)


regulator_list<-unique(df_plotreg$regulator)
regulator_list<-sort(regulator_list)
m_pi<-matrix(rep(0,length(regulator_list)*20),ncol=20)
for (i in 1:length(regulator_list)){
  x<-c(df_plotreg[df_plotreg$regulator ==regulator_list[i],3:4])
  m_pi[i,]<-c(rbind(x[[1]],x[[2]]))
}
rownames(m_pi)<-regulator_list
colnames(m_pi)<-paste0('cluster',rep(1:n_clusters,each=2),'_',rep(c('percentage','intensity'),n_clusters))
write.table(m_pi,paste0(filedir,'/regulators/regulator_sumamry.csv'),append = FALSE, sep = ",",quote=FALSE,col.names = TRUE, row.names = TRUE)

# df_sigma
num_sigma<-length(unique(df_sigma$regulator))
num_condition<-length(unique(sample_list$ConditionID))
m_sigmacondition<-matrix(0, nrow=num_sigma,ncol=num_condition,
                         dimnames = list(sort(unique(df_sigma$regulator)),unique(sample_list$ConditionID)))
m_express_condition<-matrix(0, nrow=length(gene_names),ncol=num_condition,
                            dimnames = list(gene_names,unique(sample_list$ConditionID)))
for (i in 1:num_condition){
  cur_conexp<-m_express[,which(sample_list$ConditionID==colnames(m_express_condition)[i])]
  if(is.vector(cur_conexp))
    m_express_condition[,i]<-cur_conexp
  else
    m_express_condition[,i]<-rowMeans(cur_conexp)
}
for (i in 1:num_sigma){
  sigmagenes<-df_sigma$gene[df_sigma$regulator==rownames(m_sigmacondition)[i]]
  m_sigmacondition[i,]<-colMeans(m_express_condition[sigmagenes,])

}



## rgife interface
cat('writing the clusters in transcriptional landscape into arff data file for running RGIFE\n')
df_data<-as.data.frame(t(m_express_processed))
df_data$class_label<-factor(cluster_identity)
foreign::write.arff(df_data,file='results/TLclass.arff',relation='class_label')
