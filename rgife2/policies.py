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

import sys
import itertools
import numpy as np

def get_attributes(path):
	selected_attributes = []
	for line in open (path):
		if line.startswith("RGIFE"):
			selected_size = int(line.split()[2])
			performance = float(line.split()[-1])
			continue
		if line.startswith("=="):
			continue
		selected_attributes.append(line.split()[0])
	return performance, selected_size, selected_attributes





path_results = str(sys.argv[1])
runs = int(sys.argv[2])
max_run = 0
min_run = 0
min_best = 0  # best performance for min set
max_best = 0  # best performance for max set
min_atts = 10000000000
max_atts = 0

selected_sizes = []
performances = []
selected_attributes_runs = []
union = set()
partial_union = set()

for run in range(0,runs):
	path = "{0}Run{1}/summary.txt".format(path_results, run)      
	performance, selected_size, selected_attributes = get_attributes(path)
	selected_attributes_runs.append(selected_attributes)
	selected_sizes.append(selected_size)
	performances.append(performance)
		
	for attribute in selected_attributes:
		union.add(attribute)

	if selected_size < min_atts or (selected_size == min_atts and performance > min_best ):
		min_atts = selected_size
		min_run = run
		min_best = performance
	
	if selected_size > max_atts or (selected_size == max_atts and performance > max_best ):
		max_atts = selected_size
		max_run = run
		max_best = performance
		
for run in range(0,runs):
    if len(selected_attributes_runs[run])<=30 and performances[run]>=np.mean(performances) :
        for attribute in selected_attributes_runs[run]:
            partial_union.add(attribute)

with open (path_results + "max_model.txt", "wt") as output:
	for attribute in selected_attributes_runs[max_run]:
		output.write("{0}\n".format(attribute))
output.close()

with open (path_results + "min_model.txt", "wt") as output:
	for attribute in selected_attributes_runs[min_run]:
		output.write("{0}\n".format(attribute))
output.close()

with open (path_results + "union_model.txt", "wt") as output:
	for attribute in union:
		output.write("{0}\n".format(attribute))
output.close()

with open (path_results + "partialunion_model.txt", "wt") as output:
	for attribute in partial_union:
		output.write("{0}\n".format(attribute))
output.close()

attributes=[]
frequencies = []
all_selected_attributes = list(itertools.chain.from_iterable(selected_attributes_runs))
for attribute in union:
    attributes.append(attribute)
    frequencies.append(all_selected_attributes.count(attribute))
with open (path_results + "frequencies.txt", "wt") as output:
    output.write("attribute frequency\n")
    for i in np.argsort(frequencies)[::-1]:
        output.write("{0} {1}\n".format(attributes[i],frequencies[i]))
output.close()


with open (path_results + "summarise.txt", "wt") as output:
    output.write("model size: min - {0}, max - {1}, partial union - {2}, union - {3}\n".format(len(selected_attributes_runs[min_run]),len(selected_attributes_runs[max_run]),len(partial_union), len(union)))
    output.write("run size performance\n")
    for run in range(0,runs):
        output.write("{0} {1} {2}\n".format(run,selected_sizes[run],performances[run]))
#output.close()


