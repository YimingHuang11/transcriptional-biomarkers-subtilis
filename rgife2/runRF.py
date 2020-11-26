#!/usr/bin/env python
import validation
import util
# supporting libraries
import os
import time
import sys
import numpy as np
import warnings
warnings.simplefilter("ignore")









path = str(sys.argv[1])
# dataset name (with or without arff extension)
Trainset_name = str(sys.argv[2])
Testset_name = str(sys.argv[3])


block_type, validation_schema, trees, depth, missing_values, categorical_attributes, white_list, cs_rf, cost, metric, repetitions, ordinal_attributes, different_folds, tolerance_samples, seed, cv_schema = util.parse_configuration(path+"configuration.conf")
gmean, fold_accuracy, fscore, auc = validation.runRF(path,Trainset_name, Testset_name,"BestIteration/selected_best_data.arff",trees,depth,cs_rf,cost)
with open(path+"metrics.txt", "wt") as output:
    output.write("gmean:{0}\n".format(gmean))
    output.write("accuracy:{0}\n".format(fold_accuracy))
    output.write("fscore:{0}\n".format(fscore))
    output.write("auc:{0}\n".format(auc))





