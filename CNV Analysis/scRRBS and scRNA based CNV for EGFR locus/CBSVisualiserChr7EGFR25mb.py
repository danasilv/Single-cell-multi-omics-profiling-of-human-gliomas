import pandas as pd
import os
import csv

'''path to CBS.py file'''
os.chdir('/Users/lkluegel/Documents/GBM')
import CBS

'''file to datafiles'''
os.chdir('/Users/lkluegel/Documents/GBM/CNV_data/EGFR')

samplesGBM = ["105A","105B","105C","105D",124,129,122,211,115,121]

'''import the cutoffs found by DNAcopy\
to generate these files see the commands in CNV_CBS_EGFR_25mb.R '''
outputDict = dict()
for i in samplesGBM:
    file = 'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH'+str(i)+'.csv'
    df = pd.read_csv(file)
    cutoffs = list(df['loc.end'])
    #add the first location to the list to make it compatible with previous scripts
    cutoffs = [1] + cutoffs
    #subtract one from each cutoff, as R doesn't use 0-based indexing
    cutoffs = [int(x-1) for x in cutoffs]
    outputDict['MGH'+str(i)] = cutoffs

with open('/Users/lkluegel/Documents/GBM/CNV_data/EGFR/CNV_CBS_Chr7_EGFR_25mb_All_GBM_cutoffs.csv', 'w', newline="") as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in outputDict.items():
       writer.writerow([key, value])
    
'''now open the datafiles'''
dfList  = list()
for i in samplesGBM:
    file = 'GBM_MGH'+str(i)+'_CNV_EGFR_25mb.csv'
    df = pd.read_csv(file,index_col=0)
    df = df.T
    df = df.reset_index()
    df = df.drop('index',axis=1)
    #now downsample the original dataframes by the same ratio as was used for DNAcopy
    reduce_by = 5
    df = df.groupby(df.index//reduce_by).mean()
    df.index = [x+1 for x in df.index]
    df = df.T
    df['last_col'] = ['MGH'+str(i) for x in df.index]
    dfList.append(df)
    
df = pd.concat(dfList)

CBS.CBSVisualizerGrouped(df,outputDict,'CNVs based on circular binary segmentation \n Chromosome 7 EGFR +/- 25mb','CNV_CBS_Chr7_EGFR_25mb.pdf','log')
