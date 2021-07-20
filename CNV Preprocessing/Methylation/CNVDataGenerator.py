#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#run this file as

#python CNVDataGenerator Directory MGH105A_cpgs_w20_s5.txt MGH105A_cpgs_w20_s5.txt GBM_MGH105A_CNV.csv
import sys
import os
import pandas as pd
import numpy as np

wd = sys.argv[1]
inputFile = sys.argv[2]
outputFile = sys.argv[3]
cancerType = sys.argv[4]

#change directory
os.chdir(str(wd))

#read the file
normDf = pd.read_csv(str(inputFile),sep='\t',index_col=0)

#normalise sitesDf and methDf by dividing the content of each row by the sum of the row
def rowNormalizer(inputDf):
    inputDf["rowSum"] = inputDf.sum(axis=1)
    newDf = inputDf.iloc[:,:-1].div(inputDf["rowSum"], axis = 0)
    return newDf

#run function on the dataframes
normDf = rowNormalizer(normDf)

#separate the dataset into tumor and normal tissue
groups = [cancerType,"normal"]
mask0 = [x.startswith(groups[0]) for x in list(normDf.index)]
mask1 = [x.startswith(groups[1]) for x in list(normDf.index)]
tumor = normDf[mask0]
normal = normDf[mask1]

#Copy number is given as:  CNV = cell_cpgs / median( normal_cpgs )
medianNormal = normal.median(axis = 0)

cnv = tumor.div(medianNormal)
cnv = cnv.apply(np.log2)
cnv = cnv.fillna(0)
cnv = cnv.replace(np.inf,0)
cnv = cnv.replace(-np.inf,0)

#save the file
cnv.to_csv(str(outputFile))


