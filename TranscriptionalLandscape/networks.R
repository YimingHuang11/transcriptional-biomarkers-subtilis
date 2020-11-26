################## regulons list from subtiwiki #######
######manully curated, name calibrated on Nicolas name, add repression or activation flag for each gene added####

# table(df_regulons$mode)
# dim(table(df_regulons$regulator))
# df_regulons$flag<-0
# df_regulons$flag[df_regulons$mode=='activation' | df_regulons$mode=='indirect positive regulation' | df_regulons$mode=='sigma factor' | df_regulons$mode=='positive regulation']<-1
# df_regulons$flag[grepl( 'repression', df_regulons$mode,fixed = TRUE)]<--1
# df_regulons$flag[grepl( 'negative', df_regulons$mode,fixed = TRUE)]<--1
# df_regulons$flag[grepl( 'attenuation', df_regulons$mode,fixed = TRUE)]<--1
# df_regulons$flag[grepl( 'anti-activation', df_regulons$mode,fixed = TRUE)]<--1
# df_regulons$flag[grepl( 'inhibition', df_regulons$mode,fixed = TRUE)]<--1
# df_regulons
write.table(df_regulons[,c(1:2,4:7)], "data/regulons.csv", append = FALSE, sep = ",",quote=FALSE, col.names = T, row.names = F)
##########
df_regulons<-read.table("data/regulons.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
length(table(df_regulons$regulator)) #207
length(table(df_regulons$gene)) #2288

# For cytoskype file form
#AR_list, regulon_list, sigmafactor_list, RegNetFull containing regulons and sigmafactors, RegNetlong containing all three
gene_AR<-read.table( "/data/AR_list.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
gene_regulon<-read.table("data/regulons_regulator.csv", sep = ",", head=TRUE) 
gene_sigma<-read.table( "data/regulons_sigmafactor.csv", sep = ",", head=TRUE)
# nodes table
regulonnodes<-unique(df_regulons$regulator) #207 
genenodes<-unique(df_regulons$gene) #2288
df_nodes<-data.frame(name=c(regulonnodes,genenodes),node_type=c(rep('regulator',length(regulonnodes)),rep('gene',length(genenodes))))
write.table(df_nodes, 'results/networks/v_list.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)

######regulatory network from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4796004/ ######
#####manually curated,  name calibrated on Nicolas name, regulators which are proteins set as gene names, remove * symbols##### 
df_RegNet_Jose<-read.csv("data/RegNet_JoseV2.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
dim(df_RegNet_Jose)[1] #4425
length(df_RegNet_Jose$BSU[startsWith(df_RegNet_Jose$BSU,'BSU_misc')]) #62 small RNA
length(df_RegNet_Jose$BSU[startsWith(df_RegNet_Jose$BSU,'BSU_rRNA')]) #30 ribosomal RNA
length(df_RegNet_Jose$BSU[startsWith(df_RegNet_Jose$BSU,'BSU_tRNA')]) #86 transfer RNA

regnet_regulator <- data.frame(locus.tag=character(),gene_name=character(),regulator=character(),mechanism=character(),sign=numeric(),stringsAsFactors=FALSE) 
for (i in 1:dim(gene_annotations)[1])
{
  if(gene_annotations$sigma_factor[i]!='')
  {
    sig_list<-unlist(strsplit(gene_annotations$sigma_factor[i],split = '[|]'))
    for(j in 1:length(sig_list)){
      temp_df<-data.frame(locus.tag=gene_annotations$locus[i],gene_name=gene_annotations$gene_name[i],regulator=sig_list[j],mechanism='sigma factor',sign=1,stringsAsFactors=FALSE)
      regnet_regulator = rbind(regnet_regulator,temp_df)
    }
  }
  if(gene_annotations$regulator[i]!='')
  {
    reg_list<-unlist(strsplit(gene_annotations$regulator[i],split = '[|]'))
    mech_list<-unlist(strsplit(gene_annotations$regulatory_mechanism[i],split = '[|]'))
    sign_list<-unlist(strsplit(gene_annotations$regulation_sign[i],split = '[|]'))
    for(j in 1:length(reg_list)){
      if(tolower(sign_list[j]=='nan'))
        sign_list[j]='0'
      temp_df<-data.frame(locus.tag=gene_annotations$locus[i],gene_name=gene_annotations$gene_name[i],regulator=reg_list[j],
                          mechanism=mech_list[j],sign=as.numeric(sign_list[j]),stringsAsFactors=FALSE)
      regnet_regulator = rbind(regnet_regulator,temp_df)
    }
  }
}

regnet_metabolite <- data.frame(locus.tag=character(),gene_name=character(),metabolite=character(),sign=numeric(),stringsAsFactors=FALSE) 
for (i in 1:dim(gene_annotations)[1])
{
  if(gene_annotations$metabolite[i]!='')
  {
    meta_list<-unlist(strsplit(gene_annotations$metabolite[i],split = '[|]'))
    sign_list<-unlist(strsplit(gene_annotations$metabolite_sign[i],split = '[|]'))
    for(j in 1:length(meta_list)){
      if(meta_list[j]==' @' | tolower(meta_list[j])=='nan')
        next
      meta_sublist<-unlist(strsplit(meta_list[j] ,split = ' '))
      sign_sublist<-rep('0',length(meta_sublist))
      if(sign_list[j]!=''){
        sign_sublist_new<-unlist(strsplit(sign_list[j] ,split = ' '))
        for(m in 1:length(sign_sublist_new))
          if(!is.na(sign_sublist[m]))
            sign_sublist[m]<-sign_sublist_new[m]
      }
      for(m in 1:length(meta_sublist)){
        if(tolower(sign_sublist[m])=='nan')
          sign_sublist[m]='0'
        temp_df<-data.frame(locus.tag=gene_annotations$locus[i],gene_name=gene_annotations$gene_name[i],metabolite=meta_sublist[m],sign=as.numeric(sign_sublist[m]),stringsAsFactors=FALSE)
        regnet_metabolite = rbind(regnet_metabolite,temp_df)
      }
    }
  }
}
write.table(regnet_regulator,'results/networks/regnet_regulator.csv', append = FALSE, sep = "\t",quote=FALSE,col.names = T, row.names = F)
write.table(regnet_metabolite,'results/networks/regnet_metabolite.csv', append = FALSE, sep = "\t",quote=FALSE,col.names = T, row.names = F)
########## combine the regulatory network and metabolism network###
x<-regnet_regulator_plus
colnames(x)<-c('target','source','sign','mechanism')
x<-x[,c(2,1,3,4)]
y<-regnet_metabolite[,2:4]
colnames(y)<-c('source','target','sign')
y$mechanism<-''
regnet_complete<-rbind(x,y)
write.table(regnet_complete,'results/networks/regnet_complete.csv', append = FALSE, sep = "\t",quote=FALSE,col.names = T, row.names = F)
# node list
nodenames<-unique(c(regnet_complete$source,regnet_complete$target))
nodeattrs<-rep('gene',length(nodenames))
nodeattrs[grepl('_', nodenames, fixed = TRUE) | grepl('-', nodenames, fixed = TRUE) | grepl('stringent', nodenames, fixed = TRUE) | grepl('ribo', nodenames, fixed = TRUE)]<-'Others' 
regnet_node<-data.frame(name=nodenames,nodetype=nodeattrs)
regnet_node$nodetype[regnet_node$name%in%regnet_metabolite$metabolite]<-'metabolite'
regnet_node$biomarker<-rep('NO',dim(regnet_node)[1])
regnet_node$biomarker[regnet_node$name %in% panel]<-'YES'                            
write.table(regnet_node,'results/networks/regnet_nodes.csv', append = FALSE, sep = "\t",quote=FALSE,col.names = T, row.names = F)




##################
regnet_regulator<-read.table("results/networks/regnet_regulator.csv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
regnet_regulator_plus<-read.table("results/networks/regnet_regulator_merged.csv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
regnet_metabolite<-read.table("results/networks/regnet_metabolite.csv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)
regnet_complete<-read.table("results/networks/regnet_complete.csv",quote="", sep="\t",header = TRUE, stringsAsFactors = FALSE)

dim(regnet_regulator) #5486 edges
length(unique(regnet_regulator$gene_name)) #2555
length(unique(regnet_regulator$regulator)) #259
intersect(regnet_regulator$regulator,regnet_regulator$gene_name) # 149
table(regnet_regulator$mechanism)
table(regnet_regulator$regulator)

dim(regnet_regulator_plus) #6704 edges
length(unique(regnet_regulator_plus$gene_name)) #2816
length(unique(regnet_regulator_plus$regulator)) #294
intersect(regnet_regulator_plus$regulator,regnet_regulator_plus$gene_name) # 169
table(regnet_regulator_plus$mechanism)
table(regnet_regulator_plus$regulator)

dim(regnet_metabolite) #2967 edges
length(unique(regnet_metabolite$metabolite)) #160
length(unique(regnet_metabolite$gene_name)) #1382

dim(regnet_complete) #9671 edges
length(unique(regnet_complete$source)) #1568 (294+1382)
length(unique(regnet_complete$target)) #2975 (2818+160)

#### create the graph of regulatory network#####
edge_list<-as.vector(t(as.matrix(regnet_complete[,1:2])))
g_regnet <- make_graph(edge_list,directed = FALSE)
set_edge_attr(g_regnet,name='mechanism',index=E(g_regnet),value=regnet_complete[,4])
set_edge_attr(g_regnet,name='sign',index=E(g_regnet),value=regnet_complete[,3])

panel<-scan('results/clusters/relavg/c10/biomarkers/panel.txt', what="", sep="\n")
panel1<-c("thrD", "yddJ", "yhzB", "ptsG", "S520", "murG", "yqkF", "phoR", "acuB", "yxbF")
#"yclM" "yddJ" "yhzB" "ptsG" "S520" "murG" "yqkF" "phoR" "acuB" "yxbF"
panel_slim<-panel[panel %in% c(regnet_complete$source,regnet_complete$target)]
# "yddJ" "ptsG" "murG" "yqkF" "phoR" "acuB"
panel_slim1<-panel1[panel1 %in% c(regnet_complete$source,regnet_complete$target)]

# only one shortest path
extended_g_v=list()
for(i in 1:(length(panel_slim)-1)){
  paths<-shortest_paths(g_regnet,from=panel_slim[i], to=panel_slim[(i+1):length(panel_slim)])
  extended_g_v<-append(extended_g_v,names(unlist(paths$vpath)))
}

#all shortest path
extended_g_v=list()
for(i in 1:(length(panel_slim)-1)){
  i=1
  paths<-all_shortest_paths(g_regulons,from=panel_slim[i], to=panel_slim[(i+1):length(panel_slim)])
  extended_g_v<-append(extended_g_v,names(unlist(paths$res)))
}

selected_nodes<-unique(c(unlist(extended_g_v),panel))
selected_metabolites<-c(selected_nodes,regnet_metabolite$metabolite[regnet_metabolite$gene_name%in%selected_nodes])
selected_regulators<-c(selected_nodes,regnet_regulator_plus$regulator[regnet_regulator_plus$gene_name%in%selected_nodes])
selected_nodes_plus<-unique(c(selected_nodes,selected_metabolites,selected_regulators))
write.table(selected_nodes_plus, 'results/networks/selectedgenesplus.txt', append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)

df_nodes_selected<-data.frame(name=selected_nodes_plus,biomarker=c(rep('NO',length(selected_nodes_plus))),stringsAsFactors = FALSE)

hub_cent<-hub_score(g_clusterregulon,weights=clusterregulon_netslim[,3])

dim(regnet_regulator_plus)

###############create the graph of ppi network######
df_ppi<-read.table("results/networks/PPINet_full.txt", sep=" ",header = TRUE, stringsAsFactors = FALSE)
#protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
# 4181 proteins, 1021786/2 edges
write.table(df_ppi, 'results/networks/PPINet_full.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)

df_ppi_mconfidence<-df_ppi[df_ppi$combined_score>400,]
length(df_ppi_mconfidence)
length(unique(df_ppi_mconfidence$protein1))
# 4174 proteins, 164048 edges
panel_ppiID<-paste0('224308.',gene_list$Locus_tag[match(panel,gene_list$Name)])
panel_ppiID_slim<-panel_ppiID[panel_ppiID %in% unique(df_ppi_mconfidence$protein1)]

edge_list<-as.vector(t(as.matrix(df_ppi_mconfidence[,1:2])))
g_ppi <- make_graph(edge_list,directed = FALSE)

# only one shortest path
extended_g_v=list()
for(i in 1:(length(panel_ppiID_slim)-1)){
  i=1
  paths<-shortest_paths(g_ppi,from=panel_ppiID_slim[i], to=panel_ppiID_slim[(i+1):length(panel_ppiID_slim)])
  extended_g_v<-append(extended_g_v,names(unlist(paths$vpath)))
}
selected_nodes<-unique(unlist(extended_g_v))

#for (i in 1:length(selected_nodes))
#  selected_nodes[i]<-unlist(strsplit(selected_nodes[i],split='[.]'))[2]
#selected_nodes<-gene_list$Name[match(selected_nodes,gene_list$Locus_tag)]
write.table(selected_nodes, 'results/networks/selectedgenes_mconfidence_ppi.txt', append = FALSE, sep = ",",quote=FALSE,col.names = F, row.names = F)
df_nodes_selected<-data.frame(name=selected_nodes,biomarker=c(rep('NO',length(selected_nodes))),stringsAsFactors = FALSE)
df_nodes_selected$biomarker[df_nodes_selected$name %in% panel_ppiID_slim]<-'YES'                            
write.table(df_nodes_selected, 'results/networks/selected_ppi_vlist.csv', append = FALSE, sep = ",",quote=FALSE,col.names = T, row.names = F)


##############
library(dplyr)
AR_groups<-AR_list %>% group_by(AR_id) %>% summarise(gene_lists = paste(gene_name,collapse=','))
fileConn<-file("Data/AR_list.txt")
writeLines(c(AR_groups$gene_lists,'\n'), fileConn)
close(fileConn)


hc<-hclust(dist(m_express), method = "complete", members = NULL)
mycl <- cutree(hc, k=18)
x<-names(mycl[mycl==1])
dist(m_express[x,])
