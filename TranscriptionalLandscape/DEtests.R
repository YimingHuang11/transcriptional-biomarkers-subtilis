#########################find differentially expressed genes for condition contrasts##########################################
###conditions contrasts including treatment conditions versus control conditions and similar conditions in different experiments

# matrix of design matrix  - num_samples*num_conditions: '1' treatment '0' control 'non' irrelevant to this condition
m_index<-read.table('config/design_matrix.csv', sep=",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
condition_contrasts<-colnames(m_index)

# run limma modified t-tests for all condition contrasts
# output matrix of log fold changes - num_genes*num_conditions, p-values - num_genes*num_conditions, regulation signs - num_genes*num_conditions: 1 upregualte, -1 downregulate, 0 non-significant
filedirec='results/DEGs/conditions'
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

#conditions causing few genes to response
condition_sign<-colSums(abs(m_signs))
condition_sign[condition_sign==0]
# genes invariant to all conditions
gene_sign<-rowSums(abs(m_signs))
gene_sign[gene_sign<100]

# heatmaps of sign p-values and log fold changes
x<-pheatmap::pheatmap(t(m_lfcs), fontsize_col = 0.5,show_colnames = FALSE,
         file = paste0(filedirec,"/0log2_fold_changes.pdf"),color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
         main = 'log fold changes in gene expressions under different conditions',width=12,height = 8)
pheatmap::pheatmap(t(m_sign_pvaules), fontsize_col = 0.5,cluster_rows=x$tree_row,cluster_cols=x$tree_col,show_colnames = FALSE,
                   file = paste0(filedirec,'/0sign_pvaules.pdf'),color = c("#313695","#74ADD1", "#FFFFFF", "#F46D43","#A50026"),
                   main = '-sign p-values for genes differentially under different conditions ',width=12,height = 8)
