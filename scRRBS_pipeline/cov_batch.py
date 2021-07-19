#!/usr/bin/env python

from __future__ import (division,
                        print_function,
                        absolute_import,
                        unicode_literals)

import sys
import os

cells=sys.argv[1:-2]
PDR_file=sys.argv[-2]
Description=sys.argv[-1]

counter=0
submission_list=[]

with open("QC_stats_"+Description+".csv", "w") as f:
  f.write("Cell,Methylation,PDR,Methylated_sites,Total_Reads,Biotype,Protocol,CpG_Sites,Mapping_Efficiency,CpG_sites_read,CH_conversion_rate,CG_conversion_rate\n")

for cell_name in cells:
  cov_file=cell_name+".R1.fastq_bismark_bt2_pe.bismark.cov.gz"
  PE_report=cell_name+".R1.fastq_bismark_bt2_PE_report.txt"
  submission_list.append([cov_file,cell_name,PE_report])

sp=" "

for each_cell in submission_list:
  print(each_cell)

  print("  ")
  print("  ")
  print("bash cov_practice.sh "+each_cell[0]+sp+each_cell[1]+sp+each_cell[2]+sp+PDR_file+sp+Description+" >> QC_stats_"+Description+".csv")
  print("   ")
  print("   ")
  os.system("bash cov_practice.sh "+each_cell[0]+sp+each_cell[1]+sp+each_cell[2]+sp+PDR_file+sp+Description+" >> QC_stats_"+Description+".csv")

  counter+=1
  print(counter)
  print(" ")
  print(" ")


  print("next sample")   

  
