import sys
import numpy as np

datas = []
nchar=""
ntax=""
taxas = []
for line in open(sys.argv[1]):

    if "site" in line:
        taxas = line.strip().split(",")[1:]
        continue
    line_arr = line.strip().split(',')
    ntax=len(line_arr[1:])
    data = line_arr[1:]
    datas.append("".join(data))

datas = zip(*datas)
nchar=len(datas[0])
    
print str(ntax),  str(nchar)
for tax, data in zip(taxas, datas):
    print str(tax) + "\t" + "".join(data)
