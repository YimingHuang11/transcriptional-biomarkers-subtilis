#!/usr/bin/env python
from sklearn.ensemble import RandomForestClassifier
import datasetParser
import numpy as np
import os
import sys
import subprocess

path = str(sys.argv[1])
trainfile = str(sys.argv[2])
testfile = str(sys.argv[3])
trees = int(sys.argv[4])
depth = int(sys.argv[5])


classifier = RandomForestClassifier(n_estimators = trees,max_depth=depth,class_weight = "balanced",n_jobs=4)
train_data,train_labels, attributes, class_mapping = datasetParser.read_arff_file(path+trainfile)
test_data,sample_names = datasetParser.read_expm(path+testfile)
            
#cost-sensitive learning
classifier = classifier.fit(train_data,train_labels)
#predict class for the samples
predicted_labels = classifier.predict(test_data)

with open (path+ "predicts_" + testfile.split('_')[1], "wt") as output:
    output.write("sample\tname\tpredicted label\n")
    for i in range(0,len(predicted_labels)):
        output.write("{0}\t{1}\n".format(sample_names[i],'type'+str(predicted_labels[i]+1)))
output.close()

