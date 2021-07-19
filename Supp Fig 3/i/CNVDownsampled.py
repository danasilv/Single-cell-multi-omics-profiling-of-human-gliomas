import sys
import os
import pandas as pd
import numpy as np

wd = sys.argv[1]
inputFile = sys.argv[2]
outputFile = sys.argv[3]
cancerType = sys.argv[4]
window = int(sys.argv[5])


os.chdir(str(wd))

#read the file
normDf = pd.read_csv(str(inputFile),sep='\t',index_col=0)

#normalise sitesDf and methDf by dividing the content of each row by the sum of the row
def rowNormalizer(inputDf):
    inputDf["rowSum"] = inputDf.sum(axis=1)
    newDf = inputDf.iloc[:,:-1].div(inputDf["rowSum"], axis = 0)
    return newDf

def downSampler(dataframe, window):
    newColumns = list()
    for i in range(0,dataframe.shape[1],window):
        start = i
        startChr = int(dataframe.columns[i].split(".")[0][3:])
        end = i+window
        endChr = int(dataframe.columns[i+window].split(".")[0][3:])
        while endChr > startChr:
            end = end - 1
            endChr = int(dataframe.columns[end].split(".")[0][3:])
        if end >= dataframe.shape[1]:
            end = -1
        colName = dataframe.columns[i]+"_new"
        newColumns.append(colName)
        dataframe[colName] = dataframe.iloc[:,start:end].sum(axis=1)
    dataframe = dataframe[newColumns]
    return(dataframe)
    
#run function on the dataframes
normDf = rowNormalizer(normDf)

normDf = downSampler(normDf, window)

#separate the dataset into tumor and normal tissue
groups = [cancerType,"normal"]
mask0 = [x.startswith(groups[0]) for x in list(normDf.index)]
mask1 = [x.startswith(groups[1]) for x in list(normDf.index)]
tumor = normDf[mask0]
normal = normDf[mask1]

#Copy number is given as:  CNV = cell_cpgs / mean( normal_cpgs )
#calculate the mean normal cpgs
meanNormal = normal.mean(axis = 0)

cnv = tumor.div(meanNormal)
cnv = cnv.apply(np.log2)
cnv = cnv.fillna(0)
cnv = cnv.replace(np.inf,0)
cnv = cnv.replace(-np.inf,0)

#save the file
cnv.to_csv(str(outputFile))


