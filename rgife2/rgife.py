#!/usr/bin/env python

# Copyright (C) 2015-2016 Nicola Lazzarini
# School of Computing Science, Newcastle University

# This file is part of the RGIFE heuristic.

# RGIFE is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# RGIFE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details:
# http://www.gnu.org/licenses/gpl.html


# rgife libraries
import datasetParser
import validation
import remove
import preprocessing
import util
# supporting libraries
import os
import time
import sys
import subprocess
import random
import numpy as np
import warnings
from collections import defaultdict
from sklearn.metrics import f1_score
from sklearn.ensemble import RandomForestClassifier
from operator import itemgetter
from sklearn import metrics
from sklearn.utils import resample
from sklearn.model_selection import StratifiedKFold
import pandas as pd


warnings.simplefilter("ignore")


def set_random_seeds(seed):
    if seed == -1:
        seed = int(random.random() * 100000000)
    random.seed(seed)
    np.random.seed(seed)
    validation.set_seed(seed)
    remove.set_seed(seed)
    util.set_seed(seed)
    preprocessing.set_seed(seed)
    return seed


def generate_folds(iteration, current_attributes):
    # used when multiple repetitions are required

    # current attributes is used to store inforamtion (if pre_processing>0)

    for fold in range(0, n_folds):
        subprocess.call("rm " + path_results+"it" +
                        str(actual_iteration)+"/TrainFold"+str(fold), shell=True)
        subprocess.call("rm " + path_results+"it" +
                        str(actual_iteration)+"/TestFold"+str(fold), shell=True)
    if (loocv):
        validation.loocv_folds(path_results+"it"+str(actual_iteration)+"/"+name+"_it"+str(
            actual_iteration)+".arff", path_results+"it"+str(actual_iteration))
    else:
        if pre_processing > 0:
            # re-store all the information about the CURRENT AVAILABLE ATTRIBUTES (different than the original/starting set!)
            classes[:] = []
            data_matrix[:] = []
            data_labels[:] = []
            attributes_type.clear()
            attribute_indexes.clear()
            binarised_attribute_mapping.clear()
            isCategorical.clear()
            attribute_definitions[:] = []
            relation_name, class_name = datasetParser.get_attributes_info(path_results+"it"+str(actual_iteration)+"/"+name+"_it"+str(
                actual_iteration)+".arff", attribute_definitions, attributes_type, isCategorical, attribute_indexes)

        validation.cvs_folds(path_results+"it"+str(actual_iteration)+"/"+name+"_it"+str(actual_iteration)+".arff",
                             path_results+"it"+str(actual_iteration), categorical_attributes, binarised_attribute_mapping, n_folds, db_scv)
        if pre_processing > 0:
            # pre-process the folds if necessary
            preprocessing.preprocess_dataset(pre_processing, path_results, actual_iteration,
                                             dataset_name, n_folds, binarised_attribute_mapping, isCategorical, relation_name)
            if iteration != 0:
                validation.fix_new_folds(
                    iteration, current_attributes, path_results, n_folds)
    return

def run_iteration_cv(iteration):  # accuracy, g-mean, f1-score, confusion matrix, auc, specificity and sensitivity

    rankings = {}
    iteration_metrics = {}

    current_attributes = []

    flag=0
    total_predicted_labels = []
    total_predicted_scores = []
    total_test_labels =[]

    attributes_history["it"+str(iteration)]=datasetParser.get_attributes_name(path_results+"it"+str(iteration)+"/"+name+'_it'+str(iteration)+".arff")

    for rep in range(0, repetitions):
        # probabiliites to belong to different classes for all the samples
        if different_folds == "yes":
            if rep == 0 and pre_processing > 0:
                # memorise the current list of attributes (already binarised)
                for line in open(path_results+"it"+str(iteration)+"/TrainFold0"):
                    if line.startswith("@att"):
                        current_attributes.append(line.split()[1])
            generate_folds(iteration, current_attributes)

        for fold in range(0, n_folds):

            train_data, train_labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+"it"+str(iteration)+"/TrainFold"+str(fold))
            test_data, test_labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+"it"+str(iteration)+"/TestFold"+str(fold))

            # RGIFE classifier
            # NOTE: Modify here to use a different base classifier
            classifier = RandomForestClassifier(n_estimators=int(trees), max_depth=depth, n_jobs=4)
            # cost-sensitive learning
            classifier = classifier.fit(train_data, train_labels,sample_weight=np.array([cost[i] for i in train_labels], np.float))
            # predict class for the samples
            predicted_labels = classifier.predict(test_data)
            # predict probabilities to belong to each class ((n_samples, n_classes))
            predicted_scores = classifier.predict_proba(test_data)

            if flag==0:
                total_predicted_scores = predicted_scores
                flag=1
            else:
                total_predicted_scores=np.concatenate((total_predicted_scores, predicted_scores))
            total_predicted_labels += predicted_labels.tolist()
            total_test_labels += test_labels

    ## feature importance
        for i, value in enumerate(classifier.feature_importances_):
            if attributes[i] not in rankings:
                rankings[attributes[i]] = value
            else:
                rankings[attributes[i]] += value



    sample_weight = np.array([cost[i] for i in total_test_labels], np.float)

    ## confusion matirx
    cm = metrics.confusion_matrix(total_test_labels, total_predicted_labels)
    ## accuracy (weighted average on multi classes)
    accuracy = metrics.accuracy_score(total_test_labels, total_predicted_labels, sample_weight=sample_weight)
    ## recall (weighted average on multi classes)
    recall = metrics.recall_score(total_test_labels, total_predicted_labels, average='weighted',sample_weight=sample_weight)
    ## precision
    precision = metrics.precision_score(total_test_labels, total_predicted_labels, average='weighted',sample_weight=sample_weight)
    ## fscore(weighted average on multi classes)
    fscore = metrics.f1_score(total_test_labels, total_predicted_labels, average='weighted', sample_weight=sample_weight)
    ## auc(weighted average on multi classes)
    if not loocv:
        auc_multiclass=[]
        for classno in range(0,n_classes):
            fpr, tpr, thresholds = metrics.roc_curve(total_test_labels, np.transpose(total_predicted_scores)[classno], pos_label=classno)
            auc_multiclass.append(metrics.auc(fpr, tpr))
        auc=np.average(auc_multiclass, weights=cost)


    cms.append(cm)
    # print "Confusion Matrix  Row real class, Column predicted class (combine repetitions)"
    for i in range(0, cm.shape[0]):
        row = ""
        for j in range(0, cm.shape[1]):
            row = row + "{0:.2f} ".format(cm[i][j] / repetitions)
        print (row)

    # create the file with feature importance
    with open(path_results+"it"+str(iteration)+"/attribute_score.txt", "wt") as output:
        for key, v in sorted(rankings.items(), key=itemgetter(1)):
            output.write("{0} {1}\n".format(key, rankings[key]))

    # round the performances to the 7th decimal (to fix the loss of precision problem)
    iteration_metrics["accuracy"] = round(accuracy, 7)
    iteration_metrics["fscore"] = round(fscore, 7)
    if not loocv:
        iteration_metrics["auc"] = round(auc, 7)
    iteration_metrics["recall"] = round(recall, 7)
    iteration_metrics["precision"] = round(precision, 7)

    return iteration_metrics



def run_iteration(iteration):
    iteration_metrics = run_iteration_cv(iteration)
    print ("== Metrics ==")
    for key in iteration_metrics:
        print ("{0} of iteration {1} is {2}".format(
            key, iteration, iteration_metrics[key]))
    print ("=============")
    accuracies.append(iteration_metrics['accuracy'])
    fscores.append(iteration_metrics['fscore'])
    aucs.append(iteration_metrics['auc'])
    recalls.append(iteration_metrics['recall'])
    precisions.append(iteration_metrics['precision'])
    return iteration_metrics[metric]




def all_tested(starting_index, reference_iteration, white_list,path_results):
    # check if all the attributes have been tested (removed)
    ranked = []
    for line in open(path_results+"it"+str(reference_iteration)+"/attribute_score.txt"):
        attribute = line.split()[0]
        if attribute in white_list:
            continue
        ranked.append(attribute)
    if starting_index >= len(ranked):
        return True
    else:
        return False


def check_fails(iteration, best_iteration, reference_iteration):
    # check if there exist soft fails
    best = performancies[best_iteration]
    print ("Checking previous iterations . . .")
    for i in range(reference_iteration+1, iteration+1):
        if round((best - performancies[i]), 7) <= round(tolerance, 7):
            return i
    return -1

def run_validation(best_iteration, repeats=1,alpha=10):
    current_attributes=attributes_history['it'+str(best_iteration)]

    # get reduced data
    indexs = [attributes.index(attribute) for attribute in current_attributes]
    current_data = np.asarray(data)[:,indexs]

    repetitions_accuracy=[]
    repetitions_precision=[]
    repetitions_fscore=[]
    repetitions_auc=[]
    repetitions_recall=[]
    cms=np.zeros((n_classes,n_classes))
    # repeatly generate the bootstrap resampled training set and test set
    for repeat in range(0,repeats):
        flag=0
        total_predicted_labels = []
        total_predicted_scores = []
        total_test_labels =[]
        # booststrap resampling
        data_resam,labels_resam=resample(current_data,labels,random_state=0,stratify=labels)
        skf=StratifiedKFold(n_splits=n_folds, random_state=None, shuffle=False)
        for train_index,test_index in skf.split(data_resam,labels_resam):
            train_data = np.asarray(data)[train_index,:]
            test_data = np.asarray(data)[test_index,:]
            train_labels=np.asarray(labels)[train_index]
            test_labels=np.asarray(labels)[test_index]
            classifier = RandomForestClassifier(n_estimators=int(trees), max_depth=depth, n_jobs=4)
            classifier = classifier.fit(train_data, train_labels, sample_weight=np.array([cost[i] for i in train_labels], np.float))
            predicted_labels = classifier.predict(test_data)
            predicted_scores = classifier.predict_proba(test_data)

            if flag==0:
                total_predicted_scores = predicted_scores
                flag=1
            else:
                total_predicted_scores=np.concatenate((total_predicted_scores, predicted_scores))
            total_predicted_labels += predicted_labels.tolist()
            total_test_labels += test_labels.tolist()
        sample_weight = np.array([cost[i] for i in total_test_labels], np.float)
        cm = metrics.confusion_matrix(total_test_labels, total_predicted_labels)
        accuracy = metrics.accuracy_score(total_test_labels, total_predicted_labels, sample_weight=sample_weight)
        recall = metrics.recall_score(total_test_labels, total_predicted_labels, average='weighted',sample_weight=sample_weight)
        precision = metrics.precision_score(total_test_labels, total_predicted_labels, average='weighted',sample_weight=sample_weight)
        fscore = metrics.f1_score(total_test_labels, total_predicted_labels, average='weighted', sample_weight=sample_weight)
        if not loocv:
            auc_multiclass=[]
            for classno in range(0,n_classes):
                fpr, tpr, thresholds = metrics.roc_curve(total_test_labels, np.transpose(total_predicted_scores)[classno], pos_label=classno)
                auc_multiclass.append(metrics.auc(fpr, tpr))
            auc=np.average(np.array(auc_multiclass), weights=cost)

    repetitions_accuracy.append(accuracy)
    repetitions_recall.append(recall)
    repetitions_precision.append(precision)
    repetitions_fscore.append(fscore)
    repetitions_auc.append(auc)
    cms+=cm
    #write the metrics into file
    df_metrics=pd.DataFrame(data={'Repeatition':list(range(repeats+1)),'acuracy':repetitions_accuracy,'recall':repetitions_recall,'precision':repetitions_precision,'f1_score':repetitions_fscore,'auc':repetitions_auc})
    df_metrics.to_csv(path_results+'iteration_'+str(best_iteration)+'_log.csv')

    cms/=repeats
    #print "Confusion Matrix  Row real class, Column predicted class (avg. across repetitions)"
    print(str(repeats)+"repeatitions of Bootstrap "+str(n_folds)+"-folds cross validation on the best iteration:")
    print("Average confusion Matrix")
    for i in range(0, cms.shape[0]):
        row = ""
        for j in range(0, cms.shape[1]):
            row = row + "{0:.2f} ".format(cms[i][j] / repeats)
        print (row)
    print("Metrics confidence interval (scores: "+str(100-alpha/2)+")")
    print('accuracy ['+ "{0:.4f} ".format(np.percentile(repetitions_accuracy, alpha/2))+', '+"{0:.4f} ".format(np.percentile(repetitions_accurac, 100-alpha/2))+']')
    print('fscores ['+ "{0:.4f} ".format(np.percentile(repetitions_fscore, alpha/2))+', '+"{0:.4f} ".format(np.percentile(repetitions_fscores, 100-alpha/2))+']')
    print('auc ['+ "{0:.4f} ".format(np.percentile(repetitions_auc, alpha/2))+', '+"{0:.4f} ".format(np.percentile(repetitions_auc, 100-alpha/2))+']')
    print('recall ['+ "{0:.4f} ".format(np.percentile(repetitions_recall, alpha/2))+', '+"{0:.4f} ".format(np.percentile(repetitions_recall, 100-alpha/2))+']')
    print('precision ['+ "{0:.4f} ".format(np.percentile(repetitions_precision, alpha/2))+', '+"{0:.4f} ".format(np.percentile(repetitions_precision, 100-alpha/2))+']')

    return

### read arguments from command line ####
# get the directory where configuration file and data file lie
path_results= os.getcwd()+"/"+ str(sys.argv[1])
# parse configuration file
configuration_file = str(sys.argv[2])
# dataset name (with or without arff extension)
dataset_name = str(sys.argv[3])
# name without arff extension
name = os.path.splitext(dataset_name)[0]
### set random seed ###
current_seed = set_random_seeds(seed)

### get data information ###
# number of samples
n_samples = datasetParser.get_num_samples(path_results+dataset_name)
# number of classes
n_classes = datasetParser.get_num_classes(path_results+dataset_name)
# get number of attributes
n_atts = len(datasetParser.get_attributes_name(path_results+dataset_name))
# read data
data, labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+dataset_name)


### parse configuration file ###
block_type, validation_schema, trees, depth, missing_values, categorical_attributes, white_list, cs_rf, cost, metric, repetitions, ordinal_attributes, different_folds, tolerance_samples, seed, cv_schema = util.parse_configuration(
    path_results+configuration_file)
binarised_attribute_mapping = {}

if block_type.upper() == "RBS":
    rbs = True
else:
    rbs = False

if validation_schema.upper() == "LOOCV":
    loocv = True
    n_folds = n_samples
else:
    loocv = False
    if validation_schema.upper() == "6CV":
        n_folds = 6
    else:
        n_folds = 10

if cv_schema == "db_scv":
    db_scv = True
else:
    db_scv = False

# set tolerance for good failure
# if it's a floating value then use directly the value, else value/samples
if tolerance_samples < 1:
    tolerance = tolerance_samples
else:
    tolerance = round(tolerance_samples/n_samples, 7)

if missing_values.lower() == "no" and categorical_attributes.lower() == "no":
    pre_processing = 0
if missing_values.lower() == "yes" and categorical_attributes.lower() == "yes":
    pre_processing = 1
if missing_values.lower() == "yes" and categorical_attributes.lower() == "no":
    pre_processing = 2
if missing_values.lower() == "no" and categorical_attributes.lower() == "yes":
    pre_processing = 3



# print information about the configuration and data
print ("Random Seed: {0}".format(current_seed))
print ("Configuration:")
print ("Dataset: {0}".format(name))
print ("Num Atts: {0}".format(n_atts))
print ("Num Samples: {0}".format(n_samples))
print ("Num Classes: {0}".format(n_classes))
print ("Tolerance value: {0}".format(tolerance))
print ("Missing values: {0}".format(missing_values))
print ("Categorical attributes: {0}".format(categorical_attributes))
print ("Classification cost: {0}".format(cs_rf))

if cs_rf == "auto":
    cost = n_samples/(n_classes*np.bincount(np.array(labels)))
if cs_rf == "no":
    cost = np.ones(n_classes)
for i, value in enumerate(cost):
    print ("Cost of class {0} : {1}".format(i, value))

print ("Block type: {0}".format(block_type))
print ("Validation: {0}".format(validation_schema))
print ("Repetitions of CV: {0}".format(repetitions))
print ("Different folds: {0}".format(different_folds))
print ("Performance metric: {0}".format(metric))

starting_time = time.time()

# initial block_size
block_ratio = 0.25
block_size = (int)((round)(block_ratio * n_atts))

starting_index = 0
failures = 0
attributes_history = {}
recalls = []
precisions = []
accuracies = []
fscores = []
aucs = []
cms = []
performancies = []
actual_iteration = 0
reference_iteration = 0
best_iteration = 0
previous_soft_fail = -1

subprocess.call(["mkdir", "-p", path_results + "it"+str(actual_iteration)])
subprocess.call(["cp", path_results+dataset_name, path_results+"it" +
                 str(actual_iteration)+"/"+name+"_it"+str(actual_iteration)+".arff"])


# pre-preocess the dataset
if pre_processing > 0:
    # data structures for preprocessing
    classes = []
    data_matrix = []
    data_labels = []
    attributes_type = {}
    attribute_indexes = {}

    isCategorical = {}
    attribute_definitions = []
    relation_name, class_name = datasetParser.get_attributes_info(
        dataset_name, attribute_definitions, attributes_type, isCategorical, attribute_indexes)

if (loocv):
    validation.loocv_folds(path_results+"it"+str(actual_iteration)+"/"+name+"_it"+str(
        actual_iteration)+".arff", path_results+"it"+str(actual_iteration))
else:
    validation.cvs_folds(path_results+"it"+str(actual_iteration)+"/"+name+"_it"+str(actual_iteration)+".arff",
                         path_results+"it"+str(actual_iteration), categorical_attributes, binarised_attribute_mapping, n_folds, db_scv)

if pre_processing > 0:
    preprocessing.preprocess_dataset(pre_processing, path_results, actual_iteration,
                                     dataset_name, n_folds, binarised_attribute_mapping, isCategorical, relation_name)


print ("=== Initial Iteration ===")
actual_performance = run_iteration(actual_iteration)
initial_performance = actual_performance
reference_performance = actual_performance
best_performance = actual_performance
performancies.append(actual_performance)
print ("Initial reference/best {0} is {1}".format(metric, actual_performance))

actual_iteration += 1
sys.stdout.flush()

while (True):
    start_time = time.time()
    print ("============================")
    print ("Actual iteration {0}".format(actual_iteration))
    print ("Reference iteration {0}".format(reference_iteration))
    print ("Best iteration {0}".format(best_iteration))
    print ("The best {0} is {1}".format(metric, best_performance))
    print ("The reference {0} is {1}".format(metric, reference_performance))
    print ("The block size is {0}".format(block_size))
    print ("The block ratio is {0}".format(block_ratio))
    print ("Starting index {0}".format(starting_index))

    os.system("mkdir -p " + path_results + "it"+str(actual_iteration))
    # try to remove low ranked attributes starting from starting_index
    starting_index = remove.remove_attributes(block_size, actual_iteration, reference_iteration, starting_index,
                                              white_list, path_results, name, n_folds, dataset_name, binarised_attribute_mapping, pre_processing)

    if starting_index == -1:
        print (
            "The reference dataset has not enough attributes to remove, let's reduce the block size")
        block_ratio = block_ratio * 0.25
        if (rbs):
            block_size = (int)(round(
                block_ratio * datasetParser.get_nattributes(reference_iteration, path_results, name)))
        else:
            block_size = (int)(round(block_ratio * n_atts))
        if block_size < 1:
            print ("Finish Condition")
            break
        starting_index = 0
        failures = 0
        continue

    # get the accuracy of the current iteration
    actual_performance = run_iteration(actual_iteration)
    performancies.append(actual_performance)

    if actual_performance < reference_performance:
        print ("The {0} is worse than the reference {0}!".format(metric))
        failures += 1
        print ("Consecutive Failures {0}".format(failures))
        if failures > 5:
            # check if there is a good failure, use the BEST accuracy as reference for the soft fail
            good_failure = check_fails(actual_iteration, best_iteration, max(
                previous_soft_fail, reference_iteration))

            if good_failure != -1:
                previous_soft_fail = actual_iteration
                reference_iteration = good_failure
                reference_performance = performancies[reference_iteration]
                previous_block_size = block_size
                if (rbs):
                    block_size = (int)(round(
                        block_ratio * datasetParser.get_nattributes(reference_iteration, path_results, name)))
                result = "Worse {0} but found a good failure in it {1}".format(
                    metric, good_failure)
                if block_size < 1:

                    # the block size == 1 but the current was a good iteration -> keep the block size as it was!
                    block_size = previous_block_size
                    print (
                        "Keeping the block size as it was because it was a good iteration")
                failures = 0
                starting_index = 0
                #for i in range (0, actual_iteration+1):
                 #   if (i != reference_iteration) and (i != best_iteration):
                  #      subprocess.run("rm -rf " + path_results+"it" + str(i), shell=True)
            else:
                block_ratio = block_ratio * 0.25
                if (rbs):
                    block_size = (int)(round(
                        block_ratio * datasetParser.get_nattributes(reference_iteration, path_results, name)))
                else:
                    block_size = (int)(round(block_ratio * n_atts))

                result = "Worse {0} and no good failures - > reduced block size, now {1}".format(
                    metric, block_size)
                starting_index = 0
                failures = 0
                # check stop condition
                if block_size < 1:
                    print ("Finish Condition")
                    break
        else:
            if all_tested(starting_index, reference_iteration, white_list,path_results):

                # before reducing the block_size check if there was a good failure
                good_failure = check_fails(actual_iteration, best_iteration, max(
                    previous_soft_fail, reference_iteration))

                if good_failure == -1:
                    # there was no good failure -> reduce the block size
                    block_ratio = block_ratio * 0.25
                    if (rbs):
                        block_size = (int)(round(
                            block_ratio * datasetParser.get_nattributes(reference_iteration, path_results, name)))
                    else:
                        block_size = (int)(round(block_ratio * n_atts))
                    result = "Worse {0} and no good failures -> All tested - reduced block size, now {1}".format(
                        metric, block_size)
                    # update block_size
                    starting_index = 0
                    # new one to avoid skip of block size
                    failures = 0
                    # check stop condition
                    if block_size < 1:
                        print ("Finish Condition")
                        break
                else:
                    previous_soft_fail = actual_iteration
                    reference_iteration = good_failure
                    reference_performance = performancies[reference_iteration]
                    previous_block_size = block_size
                    if (rbs):
                        block_size = (int)(round(
                            block_ratio * datasetParser.get_nattributes(reference_iteration, path_results, name)))
                    result = "All tested but found a good failure in it {0}".format(
                        good_failure)
                    if block_size < 1:
                        # the block size == 1 but the current was a good iteration -> keep the block size as it was!
                        block_size = previous_block_size
                        print (
                            "Keeping the block size as it was because it was a good iteration")
                    failures = 0
                    starting_index = 0

            else:
                result = (
                    "Worse {0}, let's remove other attributes!".format(metric))
    else:
        # check if the actual iteration is better than the BEST accuracy -> update best iteration (used for soft fails check)
        result = ""
        if actual_performance >= best_performance:
            if actual_performance == best_performance:
                result = ("The {0} is equal as the BEST! ".format(metric))
            else:
                result = ("The {0} is better than the BEST! ".format(metric))
            best_performance = actual_performance
            best_iteration = actual_iteration

        if actual_performance == reference_performance:
            result = result + \
                ("The {0} is equal as the REFERENCE!!".format(metric))
        else:
            result = result + \
                ("The {0} is better than the REFERENCE {0}!".format(metric))
        reference_iteration = actual_iteration
        reference_performance = actual_performance
        if (rbs):
            previous_block_size = block_size
            block_size = (int)(round(
                block_ratio * datasetParser.get_nattributes(reference_iteration, path_results, name)))
            if block_size < 1:
                # the block size == 1 but the current was a good iteration -> keep the block size as it was!
                block_size = previous_block_size
                print (
                    "Keeping the block size as it was because it was a good iteration")
        failures = 0
        starting_index = 0
       # for i in range (0, actual_iteration+1):
        #    if (i != reference_iteration) and (i != best_iteration):
         #       subprocess.run("rm -rf " + path_results+"it" + str(i), shell=True)
    print ("On iteration {0}: {1}".format(actual_iteration, result))
    print ("The iteration {0} took {1} minutes".format(
        actual_iteration, (time.time() - start_time)/60))
    actual_iteration += 1
    sys.stdout.flush()

print ("Summary:")
print ("The initial {0} was: {1}".format(metric, initial_performance))
print ("The best {0} was on iteration {1}: {2} with {3} attributes".format(metric, best_iteration, best_performance, len(
    datasetParser.get_attributes_name(path_results + "it"+str(best_iteration)+"/"+name + "_it"+str(best_iteration)+".arff"))))
print ("The final reference {0} was on iteration {1}: {2} with {3} attributes".format(metric, reference_iteration, reference_performance, len(
    datasetParser.get_attributes_name(path_results + "it"+str(reference_iteration)+"/"+name + "_it"+str(reference_iteration)+".arff"))))
print ("Final Best Confusion Matrix. Row real class, Column predicted class")
for i in range(0, cms[best_iteration].shape[0]):
    row = ""
    for j in range(0, cms[best_iteration].shape[1]):
        row = row + "{0:.2f} ".format(cms[best_iteration][i][j] / repetitions)
    print (row)
print ("Final Reference Confusion Matrix. Row real class, Column predicted class")
for i in range(0, cms[reference_iteration].shape[0]):
    row = ""
    for j in range(0, cms[reference_iteration].shape[1]):
        row = row + \
            "{0:.2f} ".format(cms[reference_iteration][i][j] / repetitions)
    print (row)

print ("RGIFE took {0} minutes".format((time.time() - starting_time)/60))
util.copy_final_results(path_results, best_iteration, reference_iteration,
                        pre_processing, name, binarised_attribute_mapping, dataset_name)
util.make_tarfile(path_results)
util.write_summary(path_results, metric, best_performance)
subprocess.call('rm -rf `find . -name "it*" ! -name "iterations.tar.gz"`'.format(path_results), shell=True)


attrlens_history=[len(attribute) for attribute in attributes_history.values()]
df_iterations = pd.DataFrame(data={'Iteration':list(range(actual_iteration+1)),'attr_size':attrlens_history,'acuracy':accuracies,'recall':recalls,'precision':precisions,'f1_score':fscores,'auc':aucs})
df_iterations.to_csv(path_results+'historyiterations_log.csv')

run_validation(best_iteration,repeats=50,alpha=10)
