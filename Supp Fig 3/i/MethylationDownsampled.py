import sys
import os
import pandas as pd

wd = sys.argv[1]
inputMethFile = sys.argv[2]
inputNormFile = sys.argv[3]
outputFile = sys.argv[4]
cancerType = sys.argv[5]
window = int(sys.argv[6])

#downsamples data to the desired window size given in mb
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
   
os.chdir(str(wd))

#read the file
allDf = pd.read_csv(str(inputNormFile),sep='\t',index_col=0)
methDf = pd.read_csv(str(inputMethFile),sep='\t',index_col=0)

allDf = downSampler(allDf, window)
methDf = downSampler(methDf, window)

propMethDf = methDf/allDf

#separate the dataset into tumor and normal tissue
groups = [cancerType,"normal"]
mask0 = [x.startswith(groups[0]) for x in list(propMethDf.index)]
mask1 = [x.startswith(groups[1]) for x in list(propMethDf.index)]
tumor = propMethDf[mask0]
normal = propMethDf[mask1]

meanNormal = normal.mean(axis = 0)
propTumorMethOverNormal = tumor - meanNormal
  
#save the file
propTumorMethOverNormal.to_csv(str(outputFile))
