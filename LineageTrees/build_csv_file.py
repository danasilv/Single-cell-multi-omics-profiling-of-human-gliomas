#!/usr/bin/env python
import sys
import numpy as np
import gzip as gz

"""
This file is responsible for generating the raw data file (.csv) 

0 => unmethylated
1 => methylated
- => missing

"""
filenames = sys.argv[1:]

from random import random
print >> sys.stderr, filenames, len(filenames)
sites = {}
for i, filename in enumerate(filenames):
    print >> sys.stderr, "reading", filename
    if filename[-3:] == ".gz":
        fd = gz.open(filename, 'r')
    else:
        fd = open(filename, 'r')
    for line in fd: 
        if "bedCol" in line or "start" in line: continue
        #might be where to correct for strand bias
        
        line_arr = line.strip().split()
        if "anno" in filename:
            line_arr[1] = str(int(line_arr[1]) + 1)
         
        site = line_arr[0] + '-' + line_arr[1].zfill(10)
        if site not in sites:
            sites[site] = np.asarray(["-" for _ in filenames])
    
    fd.close()

for i, filename in enumerate(filenames):
    print >> sys.stderr, "processing", filename
    if filename[-3:] == ".gz":
        fd = gz.open(filename, 'r')
    else:
        fd = open(filename, 'r')

    for line in fd:
        if "bedCol" in line or "start" in line: continue
        line_arr = line.strip().split()
        if "dan.anno" in filename:
            methylation = int(line_arr[3]) / (int(line_arr[3]) + int(line_arr[4]) +0.0)
            line_arr[1] = str(int(line_arr[1]) + 1)
        elif "bed.anno" in filename:
            methylation = line_arr[3] 
            methylation = eval (methylation.strip("'") + ".")
            line_arr[1] = str(int(line_arr[1]) + 1)
        else: 
            methylation = line_arr[3] 
            methylation = float(methylation)

        site = line_arr[0] + '-' + line_arr[1].zfill(10)
        #call missing if we have a middling methylation score
        if methylation >= .9:
            methylation = 1
        elif methylation <= .1:
            methylation = 0
        else:
            methylation = "-" 
        
        sites[site][i] = methylation
    fd.close()

print >> sys.stderr, "collecting statistics"
filenames = [filename.split("/")[-1] for filename in filenames]    

complete_data = [ len(filenames) - sites[site].tolist().count("-") for site in sites ]

stat_fd = open("stats.txt", 'w')
for i in range(1, len(filenames) + 1):
    print >> stat_fd, i, complete_data.count(i), complete_data.count(i) / (len(complete_data) + 0.0)
stat_fd.close()

print >> sys.stderr, "outputting results"

print ",".join(["site"] + filenames)

for site in sorted(sites):
    print ",".join([site] + map(str, sites[site]))


