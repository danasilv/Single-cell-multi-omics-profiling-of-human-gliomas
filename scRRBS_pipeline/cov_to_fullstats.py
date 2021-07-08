#!/usr/bin/env python

from __future__ import (division,
                        print_function,
                        absolute_import,
                        unicode_literals)

import sys
import gzip as gz
import numpy as np



if len(sys.argv) < 4 or "-h" in sys.argv or "--help" in sys.argv:
    print("cov_to_fullstats.py\t<cov.gz>\t<samplename>\t<splitting_report>\t<pdr>", file=sys.stderr)
    sys.exit()

cov = sys.argv[1]
samplename = sys.argv[2] 
report_f = sys.argv[3]
pdr_f = sys.argv[4]
btype =  sys.argv[5]
bio = sys.argv[6] 

barcodes = map(str.strip, open("barcodes.txt",'r').readlines())
barcodes = set(barcodes)

bad_bcs = open("bad_barcodes.txt",'a')
if samplename.split(".")[1] not in barcodes:
    print(samplename.split(".")[1], file=bad_bcs)
    sys.exit()
    
bad_bcs.close()
sample_meth = []


for unique_cpgs, line in enumerate(gz.open(cov, 'rt')):
    if "" in line: pass #header line
    methylation = round(float(line.strip().split()[3])/100)
    sample_meth.append(methylation)
    pass

for line in open(report_f, 'r'):
    if "Sequence pairs analysed in total" in line:
        num_reads = int(line.strip().split(":")[1].strip())
    if "Mapping efficiency" in line:
        map_effic= float(line.strip().split(":")[1].strip("%\t"))
    if "Total methylated C's in CpG context" in line: 
        meth_cpg = int(line.strip().split(":")[1].strip())
    if "Total unmethylated C's in CpG context" in line: 
        umeth_cpg = int(line.strip().split(":")[1].strip())
    if "C methylated in CHH context" in line:
        bs_rate = float(line.strip().split(":")[1].strip("%\t"))

for line in open(pdr_f, 'r'):
    if samplename in line:
        discord, cord = line.strip().split()[1:]
        PDR = float(discord) / (float(cord) + float(discord))


header="cell,methylation,PDR_total,thisMeth,total_reads,type,bio,protocol,total_cpg_no_filter,mapping_efficiency,total_cpg_gtrthan1,total_cpg_gtrthan38,avgReadCpgs_nofilter,avgReadCpgs_lessthan1CpG,avgReadCpgs_gtreql3.8CpG,bsRate,total_cpg_no_filter_down,total_reads_down,"

cols = {header: val for val, header in enumerate(header.split(','))}

output=[]

protocol= ""
filler = [""]*12
filler[-4] = bs_rate

output.append(samplename)
output.append(np.mean(sample_meth))
output.append(PDR)
output.append(sum(sample_meth))
output.append(num_reads)
output.append(btype)
output.append(bio)
output.append(unique_cpgs)
output.append(map_effic)
output.append(unique_cpgs/float(num_reads))
output.append(bs_rate)
output.append(100-float(bs_rate))


print(",".join(list(map(str, output))))
    

