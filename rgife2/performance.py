#!/usr/bin/env python


import os
import time
import sys
import numpy as np
import warnings

warnings.simplefilter("ignore")


path = str(sys.argv[1])
n_folds = int(sys.argv[2])



dic_avgmetrics =	{
  "gmean": 0.0,
  "accuracy": 0.0,
  "fscore": 0.0,
  "auc":0.0
}


for i in range(0,n_folds):
    name = path+"Fold"+str(i)+"/metrics.txt"
    dic_metrics={}
    for line in open(name):
        if line in ['\n', '\r\n', ' ']:
            continue
        key = line.split(":")[0]
        value = line.split(':')[1]
        dic_avgmetrics[key] += float(value)

        
for x,y in dic_avgmetrics.items():
    dic_avgmetrics[x]= y/n_folds

with open(path+"avg_metrics.txt", "wt") as output:
    for x,y in dic_avgmetrics.items():
        output.write("{0}:{1}\n".format(x, y))
        

    
    

