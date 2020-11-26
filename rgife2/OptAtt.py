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





# get the signatures/panel

# further manually minimising the panel
    data = False
    matrix = []
    indexes = []
    index = 0
    with open(path_results + "BestIteration/selected_best_data.arff", "wt") as output:
        output.write("@relation BestIt\n")
        for line in open(dataset_name):
            if line in ['\n', '\r\n', ' ']:
                continue
            if line.lower().startswith("@att"):
                att_name = line.split()[1]
                if att_name not in best_attributes:
                    indexes.append(index)
                    index += 1
                    continue
                index += 1
                output.write(line)
            if line.lower().startswith("@data"):
                output.write("@data\n")
                data = True
                continue
            if data:
                matrix.append(np.array(line.split(",")))
        filtered = np.delete(matrix, indexes, 1)
        for row in filtered:
            output.write(",".join(row))
#print out the current signatures used

#cross-validation
#get training data and test data

#build the RF model

#average the performance, print out





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
    
    
    
