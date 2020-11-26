#!/usr/bin/env python
import validation
import util
# supporting libraries
import os
import time
import sys
import numpy as np
import warnings
import datasetParser
warnings.simplefilter("ignore")


path = str(sys.argv[1])
n_repeats= int(sys.argv[2])
dataset_name = str(sys.argv[3])
t_frequency = int(sys.argv[4])


m_data, labels, atts, class_mapping = datasetParser.read_arff_file(path+dataset_name)
n_atts = len(atts)
m_frequencies = np.zeros(n_atts)

attributes_min=[]
attributes_max=[]
attributes_union=[]
size_max=0
size_min=10000

for i in range(0,n_repeats):
    name = path+"Rep"+str(i)+"/summary.txt"
    dic_metrics={}
    lineno=0;
    attributes_set=[]
    for line in open(name):
        if lineno>1:
            attributes_set.append(line)
            attributes_union.append(line)  
            m_frequencies[atts.index(line.rstrip())]+=1
        lineno += 1
    print("Rep {0} panel size: {1}".format(i, len(attributes_set)))
    if len(attributes_set)>size_max:
        size_max=len(attributes_set)
        attributes_max=attributes_set
    if len(attributes_set)<size_min:
        size_min=len(attributes_set)
        attributes_min=attributes_set
        


with open(path+"attributes_min.txt", "wt") as output:
    output.write("minimal {0} attributes set selected:\n".format(len(attributes_min)))   
    for x in attributes_min:
        output.write(x)
        
with open(path+"attributes_max.txt", "wt") as output:
    output.write("maximal {0} attributes set selected:\n".format(len(attributes_max)))
    for x in attributes_max:
        output.write(x)
      

attributes_union=np.unique(attributes_union)
with open(path+"attributes_union.txt", "wt") as output:
    output.write("union attributes set selected, {0} in total:\n".format(len(attributes_union)))
    for x in attributes_union:
        output.write(x)        


atts = np.array(atts)
attributes_hf = atts[np.flatnonzero(m_frequencies>=t_frequency)]
with open(path+"attributes_highfrequency.txt", "wt") as output:
    output.write("attributes selected more than {0} out of {1} times, {2} in total:\n".format(t_frequency, n_repeats, len(attributes_hf)))
    for x in attributes_hf:
        output.write(x+'\n')    
    
    

