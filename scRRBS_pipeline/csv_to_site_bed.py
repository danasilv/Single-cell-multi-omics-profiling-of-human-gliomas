#!/usr/bin/env python

from __future__ import (division,
                        print_function,
                        absolute_import,
                        unicode_literals)

import sys

if "-h" in sys.argv or "--help" in sys.argv:
    print("\nAccepts a csv and prints out a bedfile of sites. Optionally only reports sites with at least a given number of cells with data \n \n usage: python csv_to_site_bed.py [csv] [optional minimum sites] \n \n Minimum sites is minimum number of cells that must have data at any given genomic location in order to be printed into the bed. Must be an INT. Default is 1 \n\n")  
    sys.exit()

def site_generator(site_csv,min_data):
    real_data=["0","1"]
    for lines in site_csv:
        items=lines.split(",")
        if items[0]=="site":
            continue
        
        chr_loc=items[0].split("-")
        chrom=chr_loc[0]
        loc=int(chr_loc[1])
        site_data=items[1:]

        real_reads_count=0
        for each_cell in site_data:
            meth_datapoint=each_cell.strip()
            if meth_datapoint in real_data:
                real_reads_count+=1
            
        if real_reads_count>=min_data:        
            yield "\t".join(map(str, [chrom, loc, loc+1]))



csv_file=sys.argv[1]

min_site_data=1
if len(sys.argv)>2:
    min_site_data=int(sys.argv[2])

site_file=open(csv_file)
sites=site_generator(site_file, min_site_data)

for site in sites:
    print(site)


