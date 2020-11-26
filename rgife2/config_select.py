#!/usr/bin/env python
import validation
import util
# supporting libraries
import os
import time
import sys
import numpy as np
import subprocess
import warnings

warnings.simplefilter("ignore")

path_conftune = str(sys.argv[1])
path_confsets = str(sys.argv[2])
n_configs = int(sys.argv[3])
metric_name = str(sys.argv[4])



best_configID = 0
best_metric = 0.0

for i in range(0,n_configs):
#    subprocess.run(['grep', metric_name, path+"Conf"+str(i)+"/avg_metrics.txt" ])
    name = path_conftune+"Parameters_tune/Conf"+str(i)+"/avg_metrics.txt"
    dic_metrics={}
    for line in open(name):
        if line in ['\n', '\r\n', ' ']:
            continue
        key = line.split(":")[0]
        if key==metric_name:
            value = float(line.split(':')[1])
            if value > best_metric:
                best_metric = value
                best_configID = i
   
subprocess.run(['cp',path_confsets+'configuration'+str(best_configID)+'.conf',path_conftune+'configuration_best.conf'])

    
    

