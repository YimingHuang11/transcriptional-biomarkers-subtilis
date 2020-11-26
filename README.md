# transcriptional-biomarkers-subtilis

This repository contains codes used in the paper 'Computational strategies for the identification of transcriptional biomarker panels to sense cellular growth states in Bacillus subtilis'.

1. **TransciptionalLandscape**

This contains R codes used to construct transcriptional landscape which including differential expression analysis, data visualisation, data processing, transforming the transcriptomics data in UMAP, clustering of samples, and the validation of clusters and biomarkers.

Required packages: stringr, ggplot2, plotly, cowplot, RColorBrewer, plotrix, e1071, dplyr, matrixStats, data.table, igraph, pheatmap, heatmaply, ComplexHeatmap, limma, Cairo, uwot, monocle3.

please install pandoc from https://pandoc.org/installing.html if you meet error 'Saving a widget with selfcontained = TRUE requires pandoc'

Run following commands to 1) get the data pattern plots, 2) construct the Transciptional Landscape, 3) evaluate the clusters in the landscape
Rscript DataPatterns.R 
Rscript landscape.R 
Rscript LandscapeEvalulation.R <cluster identity solution file name> 


2. **RGIFE2**

This contains Python codes used to select the reduced set of genes indicative of the cluster in the transcriptional landscape. It is an update version of RGIFE http://ico2s.org/software/rgife.html. The updates include:
- Changed from Python2 to Python3
- Added the 6 folds and 8 folds cross-validation for *validation* in rgife configuration, i.e. Validation: [6CV,8CV, 10CV]
- Changed the options for *cs_rf* which decides cost sensitive learning setting for the random forest in rgife configuration as cs_rf(manual, auto, no), where 'no' means cost sensitive learning (default), 'manual' means using the user-defined class cost, and 'auto' means using the sklearn class costs calculation for auto balance: n_samples / (n_classes * np.bincount(y))
- Added the calculations for recall, precision. Made the calculation of auc suitable for multiple classes model.  
- Modified policies.py which is used to select the final solution from multiple runs of RGIFE. Apart from the min-model, max-model, union-model, we added the calculation of the optimin-model that selects the solution with smallest size and best performance and the optiunion-model that combines all the solutions smaller than threshold size. Print the summaries of multiple solutions, including the frequencies of features selected across multiple runs, the performance and panel size  of each run size.
python policies.py [path_results] [num_runs] [threshold_size]
- Added three layer nested cross-validation to report the overall performance of RGIFE heuristic model and tuning the parameters. Print the classification metrics across multiple repetitions and draw the overall roc curve for all classes.

Required libraries: NumPy, SciPy, Scikit-learn

To learn how to run and configure RGIFE2 please read the following tutorial material:
