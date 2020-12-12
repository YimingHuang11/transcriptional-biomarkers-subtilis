This repository contains codes used in the paper 'Computational strategies for the identification of transcriptional biomarker panels to sense cellular growth states in *Bacillus subtilis*', including two main parts: **Transcriptional Landscape** and **RFIFE2**.

<h1>TransciptionalLandscape</h1>

R codes used to construct transcriptional landscape which including differential expression analysis, data visualisation, data processing, transforming the transcriptomics data in UMAP, clustering of samples, and the validation of clusters and biomarkers.

Tested on R version 4.0.2 (2020-06-22) and macOS Catalina 10.15.6
Required packages: stringr, ggplot2, plotly, cowplot, RColorBrewer, plotrix, e1071, dplyr, matrixStats, data.table, igraph, pheatmap, heatmaply, ComplexHeatmap, limma, Cairo, uwot, monocle3.

please install pandoc from https://pandoc.org/installing.html if you meet error 'Saving a widget with selfcontained = TRUE requires pandoc'

To run the programme please follow the tutorial as below:

Please put the GeneExpre.csv (gene expression data with as gene names as row names, sample IDs as column names), gene_list.csv (column names: "Name", "Locus_tag", "StartV3", "EndV3", "Strand", "classif"), sample_list.csv (at least contains columns of "sampleIndex", "SampleID","ConditionID","experimentID","condition","annotation","type","medium","experiment","laboratory") under the directory of 'data', and additionally the GeneAnnotations.csv and regulons.csv (column names: "regulon", "regulator", "mode", "gene", "locus.tag", "flag")if you want to run LandscapeEvaluation.R.

In the directory of 'config' please put the design_matrix.csv for differential expression analysis, cluster_parameters.csv for scanning on clustering modle parameters and treatment_reference_ID for appointing the treatment samples and reference samples in data standardisation process.

1. **Get the data pattern plots and run differential expression tests.** In the directory of 'config' please put the design_matrix.csv for differential expression analysis. The results of log-fold-changes and p-values, are written in logfoldchanges.csv and pvalues.csv under results/DEGs, so are the differentially expressed gene lists for individual condition contrast and global heatmaps of log-fold-changes and p-values.
```
Rscript DEtests.R
```
2. **Construct the Transcriptional Landscape.** This step performs data cleaning, data standardisation, UMAP embedding learning, clustering of samples. The arguments included in the following command representing percentile threshold for high expression, percentile threshold for low expression, percentile threshold for low variance, and condition whose corresponding differentially expressed genes and related samples to be removes, are all used in the data cleaning process. In the directory of 'config' please put treatment_reference_ID.csv for appointing the treatment samples and reference samples in data standardisation process, and cluster_parameters.csv for scanning on clustering model parameters. The umap bedding result and multiple clustering solutions are written in the directory 'results/landscape'
```
Rscript landscape.R [<threshold_highp> <threshold_lowp> <threshold_varp>] [<condition_toremove>]
```
For example,
```
Rscript landscape.R
Rscript landscape.R 0.7 0.3 0.3 "sporulation late stage"
Rscript landscape.R 0.7 0.3 0.3
Rscript landscape.R "sporulation late stage"
```
3. **Evaluate the transcriptional landscape and clusters.** This step runs Differential expression for each cluster versus the others and calculates the activity levels for different sigma factors as well as transcriptional regulators in each cluster. This requires regulon information in regulons.csv to be put under 'data' directory. The arff data file required for running RGIFE is also generated here.
```
Rscript LandscapeEvalulation.R <cluster identity solution file name>
```
4. **Evaluate the biomarkers identified.** This step visualises the gene expression of biomarkers with umap embedding, heatmaps and violin plots. It also pulls out the smallest sub regulatory network that connects these biomarkers. Please put the regulatory network relevant files under the directory 'data'.
```
Rscript BiomarkerEvaluation.R <biomarker list file name> <cluster identity solution file name>
```

<h1>RGIFE2</h1>

Python codes used to select the reduced set of genes indicative of the cluster in the transcriptional landscape. It is an update version of RGIFE http://ico2s.org/software/rgife.html. The updates include:
- Changed from Python2 to Python3
- Added the 6 folds and 8 folds cross-validation for *validation* in rgife configuration, i.e. Validation: [6CV,8CV, 10CV]
- Changed the options for *cs_rf* which decides cost sensitive learning setting for the random forest in rgife configuration as cs_rf(manual, auto, no), where 'no' means cost sensitive learning (default), 'manual' means using the user-defined class cost, and 'auto' means using the sklearn class costs calculation for auto balance: n_samples / (n_classes * np.bincount(y))
- Added the calculations for recall, precision. Made the calculation of auc suitable for multiple classes model.  
- Modified policies.py which is used to select the final solution from multiple runs of RGIFE. Apart from the min-model, max-model, union-model, we added the calculation of the optimin-model that selects the solution with smallest size and best performance and the optiunion-model that combines all the solutions smaller than threshold size. Print the summaries of multiple solutions, including the frequencies of features selected across multiple runs, the performance and panel size  of each run size.

Required libraries: NumPy, SciPy, Scikit-learn, os, time, sys, subprocess, random, warnings, collections, operator, pandas, tarfile

Run RGIFE multiple times please set the working subdirectory as run$repeat_ID, for example:
```
python ./rgife.py <working directory>/run$repeat_ID/ <config file> <data file>
```
Get the statistic of multiple RGIFE runs and the solutions chosen on different policies:
```
python ./rgife/policies.py <working directory> <num_runs> <threshold_size>

```
