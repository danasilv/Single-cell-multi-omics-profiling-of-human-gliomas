#!/usr/bin/env python

from __future__ import (division,
                        print_function,
                        unicode_literals,
                        absolute_import)

import sys
from itertools import combinations
import numpy as np

header="site,methylation,PDR_total,methylation_unweighted,PDR_unweighted,thisMeth,mixedReadCount,total_reads,type,bio,protocol,total_cpg_no_filter,total_cpg_gtrthan1,total_cpg_gtrthan38,avgReadCpgs_nofilter,avgReadCpgs_lessthan1CpG,avgReadCpgs_gtreql3.8CpG,bsRate,total_cpg_no_filter_down,total_reads_down,"
cols = {header: val for val, header in enumerate(header.split(','))}

data = []
cll_data = []
normal_data = []

datasets = {}
for i, line in enumerate(open(sys.argv[1], 'r')):
    # do work
    if i == 0:
        continue
    if "methylation" in line: 
        columns = line.strip().split(',')[:11]
        columns[0] = "filename"
    else:
        row = line.strip().split(',')
        
        data.append(row)
        if int(row[7]) < 50000: #filters by min 50000
            continue   
        protocol = row[10]
        if protocol not in datasets:
            datasets[protocol]= []
        datasets[protocol].append(row)


print("cell1, cell2, average reads, average cpgs, average methylation, average conversion rate, PDR diff, protocol")

for key in datasets:
    data = np.asarray(datasets[key])
    choices = np.random.choice(range(len(data)), size=(int(len(data)/2),2), replace=False)

    for x,y in choices:
        pdr_d= abs(float(data[x][2]) - float(data[y][2]))
        total_reads_f = float(data[x][4])*.5 + float(data[y][4]) * .5
        cpgs_f = float(data[x][7])*.5 + float(data[y][7]) * .5
        methylation= float(data[x][1])*.5 + float(data[y][1]) * .5
        bs_rate_f = float(data[x][11])*.5 + float(data[y][11]) * .5
        prep = data[x][6]
        print(",".join(list(map(str, [data[x][0], data[y][0], total_reads_f, cpgs_f, methylation,bs_rate_f, pdr_d, prep]))))
        

