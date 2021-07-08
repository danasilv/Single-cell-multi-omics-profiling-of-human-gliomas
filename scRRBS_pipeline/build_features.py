#!/usr/bin/env python

from __future__ import (division,
                        print_function,
                        absolute_import,
                        unicode_literals)

import sys
import numpy as np
import gzip as gz


filenames = sys.argv[1:]

from random import random
print(filenames, len(filenames), file=sys.stderr)
sites = {}
for i, filename in enumerate(filenames):
    print("reading", filename, file=sys.stderr)
    if filename[-3:] == ".gz":
        fd = gz.open(filename, 'rt')
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
    print("processing", filename, file=sys.stderr)
    if filename[-3:] == ".gz":
        fd = gz.open(filename, 'rt')
    else:
        fd = open(filename, 'r')

    for line in fd:
        if "bedCol" in line or "start" in line: continue
        line_arr = line.strip().split()
        if "bed.anno" in filename:
            methylation = int(line_arr[3]) / (int(line_arr[3]) + int(line_arr[4]) +0.0)
            line_arr[1] = str(int(line_arr[1]) + 1)
        elif "bed2.anno" in filename:
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

print("collecting statistics", file=sys.stderr)
filenames = [filename.split("/")[-1] for filename in filenames]    

complete_data = [ len(filenames) - sites[site].tolist().count("-") for site in sites ]

stat_fd = open("stats.txt", 'w')
for i in range(1, len(filenames) + 1):
    print(i, complete_data.count(i), complete_data.count(i) / (len(complete_data) + 0.0), file=stat_fd)
stat_fd.close()

print("outputting results", file=sys.stderr)

print(",".join(["site"] + filenames))

for site in sorted(sites):
    print(",".join([site] + list(map(str, sites[site]))))


