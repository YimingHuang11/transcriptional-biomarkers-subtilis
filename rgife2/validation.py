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
import preprocessing
import datasetParser
import remove
import numpy as np
import scipy.stats as ss
import random
import sklearn.metrics.pairwise as sp
import subprocess
from sklearn.metrics import f1_score
from sklearn.ensemble import RandomForestClassifier
from operator import itemgetter
from sklearn import metrics

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    return


def cvs_folds(name, path_file, categorical_attributes, binarised_attribute_mapping, n_folds, db_scv):
    # generate the train/test folds using either: SCV or DB_SCV
    if categorical_attributes.lower() == "yes":
        samples_for_distance = []
        preprocessing.binarise_attributes(name, binarised_attribute_mapping)
        binarised_name = name.split(".")[0] + "_Binarised"
        data = False
        for line in open(binarised_name):
            if line.lower().startswith("@data"):
                data = True
                continue
            if data:
                sample = [float(i) if i !=
                          '?' else 0 for i in line.split(",")[:-1]]
                samples_for_distance.append(sample)
        samples_for_distance = np.array(samples_for_distance)
    class_mapping = {}
    class_samples = {}
    attributes = []
    samples = []
    labels = []
    data_matrix = []

    previous_line = ""
    header = ""
    data = False
    index = 0
    # read the original file and get attributes info (header to be copied in each fold)
    for line in open(name):
        if line in ['\n', '\r\n', ' ']:
            continue
        if line.startswith("%"):
            continue
        if line.lower().startswith("@rel"):
            header = header+line
            previous_line = line
            continue
        if line.lower().startswith("@att"):
            attribute = line.split()[1]
            attributes.append(attribute)
            header = header+line
            previous_line = line
            continue
        if line.lower().startswith("@data"):
            classes = [l.strip(" {}'\n\r") for l in previous_line.split(
                "{")[1].strip("}").split(",")]
            for i, label in enumerate(classes):
                class_mapping[label] = i
                class_samples[i] = []
            data = True
            continue
        if data:
            label = line.split(",")[-1].strip(" ' \n\r")
            labels.append(class_mapping[label])
            class_samples[class_mapping[label]].append(index)
            index += 1
            sample = [i if i != '?' else 0 for i in line.split(",")[:-1]]
            samples.append(sample)
            data_matrix.append(line)

    samples = np.array(samples)
    labels = np.array(samples)

    folds = {}

    for fold in range(0, n_folds):
        folds[fold] = []

    if (db_scv):
        if categorical_attributes.lower() == "no":
            samples_for_distance = samples
        neighbours = {}
        # rows that belong to a specific class
        n_samples = samples.shape[0]
        dist = sp.manhattan_distances(samples_for_distance, samples_for_distance)
        for i, sample in enumerate(samples_for_distance):
            neighbours[i] = [0]*(n_samples)
            ranked = ss.rankdata(dist[i, :], method='ordinal')
            for j, element in enumerate(ranked):
                #position = ranked[element]
                neighbours[i][int(element)-1] = j

        # db_scv
        for class_label in class_samples:
            sample_subset = class_samples[class_label]
            n = len(class_samples[class_label])//n_folds
            residuals = len(class_samples[class_label]) - (n * n_folds)
            e = random.choice(sample_subset)
            i = 0
            while (len(sample_subset) > residuals):
                folds[i].append(e)
                sample_subset.remove(e)
                # if there are no more instances
                if len(sample_subset) == 0:
                    break
                i = (i+1) % n_folds
                # get the closest samples of e
                index = 1
                closest = neighbours[e][index]
                while (closest not in sample_subset):
                    index += 1
                    closest = neighbours[e][index]
                e = closest

            # folds that will have the extra sample
            fold_extra_samples = random.sample(
                range(0, n_folds), len(sample_subset))
            for fold in fold_extra_samples:
                e = random.choice(sample_subset)
                folds[fold].append(e)
                sample_subset.remove(e)
    else:
        # scv
        for class_label in class_samples:
            # divide the sample of the class within each fold
            sample_subset = class_samples[class_label]
            n = len(class_samples[class_label])//n_folds
            for fold in range(0, n_folds):
                for i in range(0, n):
                    e = random.choice(sample_subset)
                    folds[fold].append(e)
                    sample_subset.remove(e)

            # folds that will have the extra sample
            fold_extra_samples = random.sample(
                range(0, n_folds), len(sample_subset))
            for fold in fold_extra_samples:
                e = random.choice(sample_subset)
                folds[fold].append(e)
                sample_subset.remove(e)

    # write the train/test fold files
    for fold in range(0, n_folds):
        test_indexes = folds[fold]
        with open(path_file+"/TrainFold"+str(fold), "wt") as output_train:
            with open(path_file+"/TestFold"+str(fold), "wt") as output_test:

                output_train.write(header)
                output_test.write(header)

                output_train.write("@data\n")
                output_test.write("@data\n")

                for sample in range(0, len(data_matrix)):
                    if sample in test_indexes:
                        output_test.write(data_matrix[sample])
                    else:
                        output_train.write(data_matrix[sample])


def loocv_folds(name, path_file):
    # generate the train/test folds using LOOCV
    header = ""
    data_matrix = []
    data = False

    for line in open(name):
        if line.startswith("@"):
            header = header+line
        if line.lower().startswith("@data"):
            data = True
            continue
        if data:
            data_matrix.append(line)

    n_samples = len(data_matrix)

    for fold in range(0, n_samples):
        with open(path_file + "/TrainFold"+str(fold), "wt") as output:
            output.write(header)
            for i, sample in enumerate(data_matrix):
                if i != fold:
                    output.write(sample)
        output.close()
        with open(path_file + "/TestFold"+str(fold), "wt") as output:
            output.write(header)
            output.write(data_matrix[fold])
        output.close()


def fix_new_folds(iteration, current_attributes, path_results, n_folds):
    # when having multiple REPETITONS the new folds, after being generate contains all the binarised variables from the current set of variables
    # current_attribute -> current set of attributes
    # remove all the binarised variables that are not in current_attributes
    indexes = datasetParser.get_index_atts(
        path_results + "it"+str(iteration)+"/TrainFold0")
    atts_to_remove = {}
    indx_to_remove = []
    for line in open(path_results + "it"+str(iteration)+"/TrainFold0"):
        if line.startswith("@att"):
            att = line.split()[1]
            if att not in current_attributes:
                atts_to_remove[att] = True
                indx_to_remove.append(indexes[att])
    # Remove for each fold
    for fold in range(0, n_folds):
        remove.filter_attributes_fold(path_results+"it"+str(iteration)+"/", path_results+"it"+str(
            iteration)+"/TrainFold"+str(fold), indx_to_remove, atts_to_remove, iteration)
        remove.filter_attributes_fold(path_results+"it"+str(iteration)+"/", path_results+"it"+str(
            iteration)+"/TestFold"+str(fold), indx_to_remove, atts_to_remove, iteration)
    return
    
    
def cvs_folds_outer(name, path, n_folds):  #dataname without .arff,

    class_mapping = {}
    class_samples = {}
    attributes = []
    samples = []
    labels = []
    data_matrix = []

    previous_line = ""
    header = ""
    data = False
    index = 0
    # read the original file and get attributes info (header to be copied in each fold)
    for line in open(path+name+".arff"):
        if line in ['\n', '\r\n', ' ']:
            continue
        if line.startswith("%"):
            continue
        if line.lower().startswith("@rel"):
            header = header+line
            previous_line = line
            continue
        if line.lower().startswith("@att"):
            attribute = line.split()[1]
            attributes.append(attribute)
            header = header+line
            previous_line = line
            continue
        if line.lower().startswith("@data"):
            classes = [l.strip(" {}'\n\r") for l in previous_line.split(
                "{")[1].strip("}").split(",")]
            for i, label in enumerate(classes):
                class_mapping[label] = i
                class_samples[i] = []
            data = True
            continue
        if data:
            label = line.split(",")[-1].strip(" ' \n\r")
            labels.append(class_mapping[label])
            class_samples[class_mapping[label]].append(index)
            index += 1
            sample = [i if i != '?' else 0 for i in line.split(",")[:-1]]
            samples.append(sample)
            data_matrix.append(line)

    samples = np.array(samples)
    labels = np.array(samples)

    folds = {}

    for fold in range(0, n_folds):
        folds[fold] = []
    samples_for_distance = samples
    neighbours = {}
    # rows that belong to a specific class
    n_samples = samples.shape[0]
    dist = sp.manhattan_distances(samples_for_distance, samples_for_distance)
    for i, sample in enumerate(samples_for_distance):
        neighbours[i] = [0]*(n_samples)
        ranked = ss.rankdata(dist[i, :], method='ordinal')
        for j, element in enumerate(ranked):
            #position = ranked[element]
            neighbours[i][int(element)-1] = j

    # db_scv

    for class_label in class_samples:
        sample_subset = class_samples[class_label]
        n = len(class_samples[class_label])//n_folds
        residuals = len(class_samples[class_label]) - (n * n_folds)
        e = random.choice(sample_subset)
        i = 0
        while (len(sample_subset) > residuals):
            folds[i].append(e)
            sample_subset.remove(e)
            # if there are no more instances
            if len(sample_subset) == 0:
                break
            i = (i+1) % n_folds
            # get the closest samples of e
            index = 1
            closest = neighbours[e][index]
            while (closest not in sample_subset):
                index += 1
                closest = neighbours[e][index]
            e = closest

     # folds that will have the extra sample
        fold_extra_samples = random.sample(range(0, n_folds), len(sample_subset))
        for fold in fold_extra_samples:
            e = random.choice(sample_subset)
            folds[fold].append(e)
            sample_subset.remove(e)

    # write the train/test fold files
    for fold in range(0, n_folds):
        test_indexes = folds[fold]
        subprocess.run(['mkdir','-p', path+'Fold'+str(fold)])
        with open(path+'Fold'+str(fold)+"/"+name+"_TrainFold"+".arff", "wt") as output_train:
            with open(path+"Fold"+str(fold)+"/"+name+"_TestFold"+".arff", "wt") as output_test:

                output_train.write(header)
                output_test.write(header)

                output_train.write("@data\n")
                output_test.write("@data\n")

                for sample in range(0, len(data_matrix)):
                    if sample in test_indexes:
                        output_test.write(data_matrix[sample])
                    else:
                        output_train.write(data_matrix[sample])
                        
                        
                        
#run random forest with selected attributes
def runRF_new(path_results,train_file, test_file,selected_attributes,trees,depth,cs_rf,cost): 
   
    # number of samples
    n_samples = datasetParser.get_num_samples(path_results+selectedAtt_file)
    # number of classes
    n_classes = datasetParser.get_num_classes(path_results+selectedAtt_file)
    
    # return the indexs of selected attributes    
    classifier = RandomForestClassifier(n_estimators=int(trees), max_depth=depth,class_weight = "balanced",n_jobs=4)
    train_data, train_labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+train_file)
    test_data, test_labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+test_file)
    indexs = attributes.index(selected_attributes)
    
    #get the reduced data
    attributes = np.asarray(attributes)[indexs]
    attributes.append('stress_type')
    train_data = np.asarray(train_data)[:,indexs]
    test_data = np.asarray(test_data)[:,indexs]

    # performance variables
    correct_classified = 0
    total_samples = 0
    true_labels = []

    misclassified_samples=[]
    cm = np.zeros((n_classes, n_classes))
    
    # cost-sensitive learning
    if cs_rf == "yes":
        classifier = classifier.fit(train_data, train_labels, sample_weight=np.array([cost[i] for i in train_labels], np.float))
    else:
        classifier = classifier.fit(train_data, train_labels)

    # predict class for the samples
    predicted_labels = classifier.predict(test_data)
    # predict probabilities to belong to each class
    predicted_percentages = classifier.predict_proba(test_data)
    true_labels += test_labels

    # get the probabilities for each class, calculate the auc for each class vs others and then average
    auc = 0   
    for n in range(0,n_classes):
        probabilities = []
        for i, value in enumerate(predicted_percentages):
            for k, pr in enumerate(value):
                if k==n:
                    probabilities.append(pr)  
        fpr, tpr, thresholds = metrics.roc_curve(test_labels, probabilities, pos_label=n)
        auc += metrics.auc(fpr, tpr)
    auc = auc / n_classes 
    
    # calculate gmean fold
    countsReal = {}
    countsCorr = {}
    for cl in class_mapping.values():
        countsReal[cl] = countsCorr[cl] = 0 
    for i, label in enumerate(predicted_labels):
        countsReal[test_labels[i]] += 1
        if label == test_labels[i]:
            countsCorr[test_labels[i]] += 1
        cm[test_labels[i]][label] += 1
        
        gmean = 1
        for cl in countsReal.keys():
            if countsReal[cl] > 0:
                rate = countsCorr[cl]/countsReal[cl]
                if rate == 0:
                    rate = 0.1/n_samples
            else:
                rate = 1
            gmean *= rate
        gmean = pow(gmean, 1/n_classes)

        fold_correct = 0
        for i, label in enumerate(predicted_labels):
            total_samples += 1
            if label == test_labels[i]:
                correct_classified += 1
                fold_correct += 1

        fold_accuracy = fold_correct / len(test_labels)

        # calculate fscore (get fscore for each class and then average)
        fscore = np.mean(f1_score(test_labels, predicted_labels, average=None))
        
    for i in range(0, cm.shape[0]):
        row = ""
        for j in range(0, cm.shape[1]):
            row = row + "{0:.2f} ".format(cm[i][j])
        print (row)
         
    return gmean, fold_accuracy, fscore, auc
    
def runRF(path_results,train_file, test_file,selectedAtt_file,trees,depth,cs_rf,cost): 
   
    # number of samples
    n_samples = datasetParser.get_num_samples(path_results+selectedAtt_file)
    # number of classes
    n_classes = datasetParser.get_num_classes(path_results+selectedAtt_file)
    
    # return the indexs of selected attributes    
    classifier = RandomForestClassifier(n_estimators=int(trees), max_depth=depth,class_weight = "balanced",n_jobs=4)
    train_data, train_labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+train_file)
    test_data, test_labels, attributes, class_mapping = datasetParser.read_arff_file(path_results+test_file)
    indexs = datasetParser.get_index_selectedatts(path_results+selectedAtt_file,attributes)
    
    #get the reduced data
    attributes = np.asarray(attributes)[indexs]
    indexs = indexs[:-1] #delete class label line
    train_data = np.asarray(train_data)[:,indexs]
    test_data = np.asarray(test_data)[:,indexs]

    # performance variables
    correct_classified = 0
    total_samples = 0
    true_labels = []

    misclassified_samples=[]
    cm = np.zeros((n_classes, n_classes))
    
    # cost-sensitive learning
    if cs_rf == "yes":
        classifier = classifier.fit(train_data, train_labels, sample_weight=np.array([cost[i] for i in train_labels], np.float))
    else:
        classifier = classifier.fit(train_data, train_labels)

    # predict class for the samples
    predicted_labels = classifier.predict(test_data)
    # predict probabilities to belong to each class
    predicted_percentages = classifier.predict_proba(test_data)
    true_labels += test_labels

    # get the probabilities for each class, calculate the auc for each class vs others and then average
    auc = 0   
    for n in range(0,n_classes):
        probabilities = []
        for i, value in enumerate(predicted_percentages):
            for k, pr in enumerate(value):
                if k==n:
                    probabilities.append(pr)  
        fpr, tpr, thresholds = metrics.roc_curve(test_labels, probabilities, pos_label=n)
        auc += metrics.auc(fpr, tpr)
    auc = auc / n_classes 
    
    # calculate gmean fold
    countsReal = {}
    countsCorr = {}
    for cl in class_mapping.values():
        countsReal[cl] = countsCorr[cl] = 0 
    for i, label in enumerate(predicted_labels):
        countsReal[test_labels[i]] += 1
        if label == test_labels[i]:
            countsCorr[test_labels[i]] += 1
        cm[test_labels[i]][label] += 1
        
        gmean = 1
        for cl in countsReal.keys():
            if countsReal[cl] > 0:
                rate = countsCorr[cl]/countsReal[cl]
                if rate == 0:
                    rate = 0.1/n_samples
            else:
                rate = 1
            gmean *= rate
        gmean = pow(gmean, 1/n_classes)

        fold_correct = 0
        for i, label in enumerate(predicted_labels):
            total_samples += 1
            if label == test_labels[i]:
                correct_classified += 1
                fold_correct += 1

        fold_accuracy = fold_correct / len(test_labels)

        # calculate fscore (get fscore for each class and then average)
        fscore = np.mean(f1_score(test_labels, predicted_labels, average=None))
        
    for i in range(0, cm.shape[0]):
        row = ""
        for j in range(0, cm.shape[1]):
            row = row + "{0:.2f} ".format(cm[i][j])
        print (row)
         
    return gmean, fold_accuracy, fscore, auc
