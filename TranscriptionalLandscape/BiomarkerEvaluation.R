args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Cluster identity input file name and biomarker list file name must be supplied", call.=FALSE)
}
######## read biomarkers #######
biomarkerfile<-args[1]
panel<-scan(biomarkerfile, what="", sep="\n")

##### read data ########
source("data.R")
source("util.R")
cat('Reading processed data\n')
m_express_rel<-as.matrix(read.table('data/m_express_rel.csv',sep=",", header = TRUE, check.names = FALSE,stringsAsFactors = FALSE))
genenames_reduced<-rownames(read.table('data/selectedGenes.csv',sep=",", header = TRUE, stringsAsFactors = FALSE))

index_treatment<-sample_list$sampleIndex[sample_list$SampleID %in% colnames(m_express_rel)]
m_express_processed<-m_express_rel[genenames_reduced,]

cat('Reading cluster solution\n')
clusterfile<-args[2]
cluster_identity<-read.table(clusterfile,sep=",", header = FALSE, stringsAsFactors = FALSE)
cluster<-cluster_identity[,2]
names(cluster)<-cluster_identity[,1]
n_clusters<-length(unique(cluster))

filedir<-'results/biomarkers'
if (!dir.exists(filedir)) {dir.create(filedir)}

##### smallest connected sub regulatory network with the biomarkers #####
regnet_complete<-read.table("data/RegNet.tsv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
regnet_regulator<-read.table("data/regnet_regulator.tsv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
regnet_metabolite<-read.table("data/regnet_metabolite.tsv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)

panel_slim<-panel[panel %in% c(regnet_complete$source,regnet_complete$target)]

#### create the graph of regulatory network#####
edge_list<-as.vector(t(as.matrix(regnet_complete[,1:2])))
g_regnet <- make_graph(edge_list,directed = FALSE)
set_edge_attr(g_regnet,name='mechanism',index=E(g_regnet),value=regnet_complete[,'mechanism'])
set_edge_attr(g_regnet,name='sign',index=E(g_regnet),value=regnet_complete[,'sign'])

# only one shortest path
extended_g_v=list()
for(i in 1:(length(panel_slim)-1)){
  paths<-shortest_paths(g_regnet,from=panel_slim[i], to=panel_slim[(i+1):length(panel_slim)])
  extended_g_v<-append(extended_g_v,names(unlist(paths$vpath)))
}

selected_nodes<-unique(c(unlist(extended_g_v),panel))
selected_metabolites<-c(selected_nodes,regnet_metabolite$metabolite[regnet_metabolite$gene_name%in%selected_nodes])
selected_regulators<-c(selected_nodes,regnet_regulator$regulator[regnet_regulator$gene_name%in%selected_nodes])
selected_nodes_plus<-unique(c(selected_nodes,selected_metabolites,selected_regulators))
write.table(selected_nodes_plus, 'results/biomarkers/biomarkers_supernodes.txt', append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)


##### visualisation of biomarker expressions ########
### heatmaps
expsigncluster<-function(clusterexps){
  return (sum(clusterexps>0)/length(clusterexps))
}


df_expressrel<-as.data.frame(t(m_express_rel))
df_expressrel$condition<-sample_list$ConditionID[index_treatment]
df_expressrel$cluster<-cluster
df_expressrel$conditioncl<-paste0(df_expressrel$condition,'_',df_expressrel$cluster)

# biomarkers*conditions
cat('Plot the heatmap of average biomarker expressions per condition\n')
df_exppercondition<-aggregate(df_expressrel, list(df_expressrel$condition), mean)
rownames(df_exppercondition)<- df_exppercondition$Group.1
df_conditionpercluster<-data.frame(condition=sample_list$ConditionID[match(names(cluster),sample_list$SampleID)],cluster=cluster)
df_conditionpercluster$condition<-df_expressrel$conditioncl

# average expression per cluster
cat('Plot the heatmap of average biomarker expressions per cluster\n')
df_exppercluster<-aggregate(df_expressrel[genenames_reduced], list(df_expressrel$cluster), mean)
df_expperconditioncl<-aggregate(df_expressrel[genenames_reduced], list(df_expressrel$conditioncl), mean)
df_expsigncluster<-aggregate(df_expressrel[genenames_reduced], list(df_expressrel$cluster), expsigncluster)


rownames(df_expperconditioncl)<- df_expperconditioncl$Group.1
df_expperconditioncl$cluster<-df_conditionpercluster$cluster[match(df_expperconditioncl$Group.1,df_conditionpercluster$condition)]
df_expperconditioncl$cluster<-paste0('cluster',df_expperconditioncl$cluster)
df_expperconditioncl<-df_expperconditioncl[order(df_expperconditioncl$cluster),]
anno_row<-data.frame(row.names=rownames(df_expperconditioncl),cluster=df_expperconditioncl$cluster)
pheatmap::pheatmap(df_expperconditioncl[,panel],  cluster_rows = F,show_colnames = T, show_rownames = T,
         annotation_row  =anno_row,angle_col = 0,
         file = paste0(filedir,"/heatmap_conditionpercluster.pdf"),width=8,height=16)

ht1=ComplexHeatmap::Heatmap(df_expperconditioncl[,panel],row_split=df_expperconditioncl$cluster,
                            row_names_side='left',column_dend_side='bottom',column_names_rot = 45,row_names_gp = gpar(fontsize = 7),
                            column_names_gp = gpar(fontsize = 10),width=10,height=12,name='Relative Gene Expression',
                            column_names_centered=TRUE,row_names_centered=TRUE)
ht2=draw(ht1)
rownames(df_exppercluster)<- df_exppercluster$Group.1
# column order adjust to be the same
order_clusterno<-as.numeric(substring(names(row_order(ht2)),8))
order_geneindex<-column_order(ht2)
ComplexHeatmap::Heatmap(df_expsigncluster[order_clusterno,panel[order_geneindex]],cluster_rows=FALSE,cluster_columns=FALSE,
                        row_names_side='left',column_names_side='bottom',column_names_rot = 45,row_names_rot = 90,
                        row_names_gp = gpar(fontsize = 17),column_names_gp = gpar(fontsize = 20),
                        width=10,height=15,name='Relative Gene Expression',row_labels = paste0('cluster',order_clusterno),
                        column_names_centered=TRUE,row_names_centered=TRUE)
pheatmap::pheatmap(df_expsigncluster[order_clusterno,panel[order_geneindex]], show_colnames = T, cluster_rows = F, cluster_cols = F,show_rownames = T,angle_col = '0',
         file = paste0(filedir,"/hm_clusterpercent.pdf"),width=10,height=15)

pheatmap::pheatmap(df_exppercluster[order_clusterno,panel[order_geneindex]], show_colnames = T, cluster_rows = F,cluster_cols = F,show_rownames = T,angle_col = '0',
         file = paste0(filedir,"/hm_clusterintensity.pdf"),width=8,height=8)
write.table(panel[order_geneindex],paste0(filedir,"/panel_sorted.txt"),append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)


### violin plots
filedire<-paste0(filedir,'/violinplots')
if (!dir.exists(filedire)) {dir.create(filedire)}

cat('Plot the violin plots of biomarker expressions in different clusters\n')
minexp<-min(m_express_rel)
maxexp<-max(m_express_rel)
for (i in order_clusterno){
  data<-as.data.frame(as.table(m_express_rel[panel[order_geneindex],names(which(cluster==i))]))
  colnames(data)<-c('biomarker','sample','expression')

  g<-data %>%
    ggplot(aes(x=biomarker, y=expression, fill=biomarker)) +
    geom_violin(show.legend = FALSE) +
    scale_fill_manual(values = brewer.pal(10,'Set3') )+
    theme(legend.position="none") +ylim(minexp, maxexp)+ theme_bw() +ggtitle("") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_blank(),axis.text.x=element_blank())
    ggsave(paste0(filedire,'/violin_',i,'.pdf'),g, height = 5, width = 10)
}

### the expression scale of individual biomarkers on umap
cat('Indicate the biomarker expressions by colour in orignial umap embedding\n')

umap_embedding<-read.table('results/landscape/umap_embedding.csv',sep=',',stringsAsFactors = FALSE)
m_express_umap<-as.data.frame(umap_embedding)
colnames(m_express_umap)<-c('component1','component2')
m_express_umap<-cbind(m_express_umap,sample_list[index_treatment,])
m_express_umap$cluster<-as.factor(cluster)

filedire<-paste0(filedir,'/umapbiomarkers')
if (!dir.exists(filedire)) {dir.create(filedire)}
setwd(filedire)

for (i in 1:length(panel)){
  m_express_umap[,panel[i]]<-m_express_rel[panel[i],]
  p<-ggplot(m_express_umap,aes(x=component1,y=component2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+geom_point(size=4,aes_string(colour=panel[i]))+theme(legend.title =element_text(size=30),legend.text =element_text(size=16),axis.text =element_text(size=20),axis.title =element_text(size=30))+scale_color_gradient(low="blue", high="red")
  ggsave(file=paste0("umap_",panel[i],".pdf"), plot = p,units="in", width=10, height=8, dpi = 300)
  p<-ggplotly(p)
  htmlwidgets::saveWidget(p, file=paste0("umap_",panel[i],".html"))
}
p<-ggplot(m_express_umap,aes(x=component1,y=component2,text =paste("Sample ID:", SampleID,"\nannotation:", annotation,"\nExperiment ",experiment)))+geom_point(size=4,aes(colour=cluster))+scale_color_manual(values = custom.col[1:n_clusters])+theme(legend.title =element_text(size=30),legend.text =element_text(size=16),axis.text =element_text(size=20),axis.title =element_text(size=30))
ggsave(file="umap_cluster.pdf", plot = p,units="in", width=10, height=8, dpi = 300)
p<-ggplotly(p)
htmlwidgets::saveWidget(p, file="umap_cluster.html")

setwd('../../../')
