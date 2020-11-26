#!/usr/bin/env python
import validation
# supporting libraries
import os
import subprocess
import warnings
warnings.simplefilter("ignore")
import sys

#arguments --path (e.g Exp1/), data file (e.g data.arff), n_folds1, n_folds2,n_parametersets
path = str(sys.argv[1])
dataset_name = str(sys.argv[2])
n_folds1 = int(sys.argv[3])
n_folds2 = int(sys.argv[4])
n_parametersets = int(sys.argv[5])
# name without arff extension
name = os.path.splitext(dataset_name)[0]


#To test the performance of this model can achieve
# First layer: K1-folds cross validtion for report the final performance 
for i in range(0,n_folds1):  
    validation.cvs_folds_outer(name, path, n_folds1)       # training data folds store as path+'Fold/'+str(fold)+"/"+name+"_TrainFold.arff"
    # Second layer: n configuration setttings
    for j in range (0,n_parametersets):
        subprocess.run(['mkdir','-p', path+'Fold'+str(i)+'/Parameters_tune/Conf'+str(j)])
        subprocess.run(["cp", path+'Fold'+str(i)+"/"+name+"_TrainFold"+".arff", path+'Fold'+str(i)+"/Parameters_tune/Conf"+str(j)+"/"+dataset_name])
        validation.cvs_folds_outer(name, path+'Fold'+str(i)+"/Parameters_tune/Conf"+str(j)+"/", n_folds2)
        subprocess.run(["rm",  path+'Fold'+str(i)+"/Parameters_tune/Conf"+str(j)+"/"+dataset_name])
        # Third layer: K2-folds cross validation to report the performance under certain configuration
        for m in range (0, n_folds2):
            subprocess.run(["cp", path+"config_sets/configuration"+str(j)+".conf", path+'Fold'+str(i)+"/Parameters_tune/Conf"+str(j)+"/Fold"+str(m)+"/configuration.conf"])  #copy the configuration file to each config_j/fold_m

# To choose the best configuration for running the model on whole dataset
# First layer: n configuration setttings             
for j in range (0,n_parametersets):            
    # Second layer: K-fold cross validation
    for i in range(0,n_folds1):
        subprocess.run(['mkdir','-p', path+'Parameters_tune/Conf'+str(j)+'/Fold'+str(i)])
        subprocess.run(["cp", path+"config_sets/configuration"+str(j)+".conf", path+"Parameters_tune/Conf"+str(j)+"/Fold"+str(i)+"/configuration.conf"]) #copy the configuration file to each Parameters_tune/config_j/fold_i  
        subprocess.run(['cp', path+"Fold"+str(i)+"/"+name+"_TrainFold.arff", path+"Parameters_tune/Conf"+str(j)+"/Fold"+str(i)+"/"+name+"_TrainFold.arff"])  #copy the data
        subprocess.run(['cp', path+"Fold"+str(i)+"/"+name+"_TestFold.arff", path+"Parameters_tune/Conf"+str(j)+"/Fold"+str(i)+"/"+name+"_TestFold.arff"])                   

