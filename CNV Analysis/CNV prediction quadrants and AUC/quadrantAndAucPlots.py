import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve
from statsmodels.stats.inter_rater import cohens_kappa 
from mpl_toolkits.mplot3d import Axes3D
import math
from sklearn.utils import resample
from scipy.stats import ttest_ind
import seaborn as sns
import time
import os
import csv
from scipy.interpolate import InterpolatedUnivariateSpline 
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy import stats
import math

#downsamples methylation data for different window sizes
def downSampleRawFile(df,resample):
    cols = df.columns.values[::resample]
    for i in range(len(cols)):
        if i < len(cols)-1:
            df[cols[i]] = df.loc[:,cols[i]:df.columns.values[list(df.columns.values).index(cols[i])+resample-1]].sum(axis=1)
        else:
            df[cols[i]] = df.loc[:,cols[i]:].sum(axis=1)
    df = df[cols]
    return(df)
    
#find the column names of the start and ends of chromosomes
def getChromosomeCutoffs(df):
    starts = list()
    ends = list()
    columnNames = df.columns
    a = [x.split(".") for x in columnNames]
    for i in range(len(a)-1):
        rightChr = int(a[i+1][0][0:])
        leftChr = int(a[i][0][0:])
        if leftChr < rightChr:
            starts.append(columnNames[i+1])
            ends.append(columnNames[i])
    starts.insert(0,columnNames[0])
    ends.append(columnNames[-1])
    return(starts,ends)
   
#row normalize the methylation data
def rowNormalizer(inputDf):
    inputDf["rowSum"] = inputDf.sum(axis=1)
    newDf = inputDf.iloc[:,:-1].div(inputDf["rowSum"], axis = 0)
    return(newDf)
  
#get methylation CNV by chromosome
def rawMethFileConverterChr(wd,inputFile,outputFile,cancerType):
    os.chdir(str(wd))
    
    #read the file
    normDf = pd.read_csv(str(inputFile),sep='\t',index_col=0)    
    normDf.columns = [x[3:] for x in normDf.columns.values]

    starts,ends = getChromosomeCutoffs(normDf)
    chromosomes = list(set([int(str(x).split('.')[0]) for x in normDf.columns.values]))
    chromosomes.sort()
    
    for i in range(len(chromosomes)):
        start = starts[i]
        end = ends[i]
        start_index = list(normDf.columns.values).index(start)
        end_index = list(normDf.columns.values).index(end)
        normDf[chromosomes[i]] = normDf.iloc[:,start_index + 1:end_index - 1].sum(axis=1)

    normDf = normDf.iloc[:,-len(chromosomes):]
    
    normDf.columns = [str(x)+'.'+str(x) for x in normDf.columns.values]
    #normalise sitesDf and methDf by dividing the content of each row by the sum of the row
    normDf = rowNormalizer(normDf)
    
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

#import methylation files
samplesGBM = ["105A","105B","105C",124,122]

for sample in samplesGBM:
    rawMethFileConverterChr('w1_s1','MGH%s_cpgs_w1_s1.txt' %sample,'MGH%s_cpgs_whole_chromosomes.csv' %sample,'GBM')

#convert raw scRRBS to CNV for different windows
def rawMethFileConverter(wd,inputFile,outputFile,cancerType,resample):
    os.chdir(str(wd))
    
    #read the file
    normDf = pd.read_csv(str(inputFile),sep='\t',index_col=0)
    normDf.columns = [x[3:] for x in normDf.columns.values]
    #get chromosome cutoffs
    starts,ends = getChromosomeCutoffs(normDf)
    dfList = list()
    for i in range(len(starts)):
        try:
            curDf = normDf.loc[:,starts[i]:normDf.columns.values[list(normDf.columns.values).index(ends[i])+1]]
        except:
            curDf = normDf.loc[:,starts[i]:]
        curDf = downSampleRawFile(curDf,resample)
        dfList.append(curDf)
    
    normDf = pd.concat(dfList,axis=1)
 
    #normalise sitesDf and methDf by dividing the content of each row by the sum of the row
    normDf = rowNormalizer(normDf)
    
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

#generate the resampled methylation CNV data
#the scRRBS data is using a 1mb window with 1mb slide
resamples = [1,2,3,5,10,20,30,40,50]

for z in resamples:
    for sample in samplesGBM:
        rawMethFileConverter('w1_s1','MGH%s_cpgs_w1_s1.txt' %sample,'MGH%s_cpgs_w%s_s%s.csv' %(sample,z,z),'GBM',z)

#do the same for more granular data (using 0.1mb windows and slide)
resamplesFiner = [1,2,3,4,5]

for z in resamplesFiner:
    for sample in samplesGBM:
        rawMethFileConverter('w0.1_s0.1','MGH%s_cpgs_w0.1_s0.1.txt' %sample,'MGH%s_cpgs_w%s_s%s.csv' %(sample,'0.'+str(z),'0.'+str(z)),'GBM',z)


#import the scRRBS CNV files for each sample in the list
def methCNVFileImporter(sampleList,fileNameEnd):
    methylationCNVFiles = ["MGH"+str(x)+fileNameEnd for x in samplesGBM]
    methCNVDfs = [pd.read_csv(x,index_col=0) for x in methylationCNVFiles]
    methCNVDf = pd.concat(methCNVDfs)
    return(methCNVDf)

#the windows are now out of order (alphabetical order), so split the index and resort them
def process_index(k):
    return (int(k.split(".")[0]),int(k.split(".")[1]))

#run the scripts for scRRBS 
os.chdir('w1_s1')

for z in resamples:
    methCNVDf = methCNVFileImporter(samplesGBM,'_cpgs_w%s_s%s.csv' %(z,z))

    #filter out cells that did not pass QC
    methCNVFilterFile = 'allQC_DNAme_POST_filtered_50K_09112019.csv'
    methCNVFilterDf = pd.read_csv(methCNVFilterFile,index_col=0)
    methCNVFilter = list(methCNVFilterDf['Cell'])
    for i in range(len(methCNVFilter)):
        if methCNVFilter[i][:4] != 'GBM_':
            methCNVFilter[i] = 'GBM_' + methCNVFilter[i]
        
    for row in methCNVDf.iterrows():
        if '.'.join(row[0].split('.')[:2]) not in methCNVFilter:
            methCNVDf = methCNVDf.drop(row[0])
    
    methCNVDf.columns = [methCNVDf.columns.values[x].split('.')[0] + '.' + str(x+1) for x in range(len(methCNVDf.columns.values))]
        
    #import the scRNA CNV data
    samplesScRNA = ["105A","105B","105C",124,122]
    scRNADfList = list()
    
    #import the scRNA based CNV data
    for i in samplesScRNA:
        file = 'inferCNV_MGH%s/preliminary_modified_expression.csv' %i
        df = pd.read_csv(file,index_col=0)
        df = df.T
        annotationFile = 'inferCNV_MGH%s/MGH%s_annotation_file_inferCNV.txt' %(i,i)
        annoDf = pd.read_csv(annotationFile,sep='\t',index_col=0,header=None)
        df = df.join(annoDf)
        df = df[df.iloc[:,-1] == 'tumor']
        scRNADfList.append(df)
    
    #get the common genes in each sample
    common_cols = list(set.intersection(*(set(df.columns) for df in scRNADfList)))
    #concatenate using only the common columns
    scRNACNVDf = pd.concat([df[common_cols] for df in scRNADfList])
    scRNACNVDf = scRNACNVDf.drop(1,axis=1)
    
    #filter out cells that did not pass QC
    scRNAQCFilterFile = 'allQC_POST_filtered_scRNA_09122019.csv'
    scRNAQCFilterDf = pd.read_csv(scRNAQCFilterFile,index_col=0)
    scRNAQCFilter = list(scRNAQCFilterDf['X1'])
    for row in scRNACNVDf.iterrows():
        if row[0] not in scRNAQCFilter:
            scRNACNVDf = scRNACNVDf.drop(row[0])
            
    #find the window start and ends for all windows in the methylation file
    startsNew,endsNew = getChromosomeCutoffs(methCNVDf)
    windowWidth = z
    windowStarts = list()
    chromosomes = [int(x.split('.')[0][0:]) for x in methCNVDf.columns.values]
    
    for col in methCNVDf.columns.values:
        if col in startsNew:
            windowStart = 0
        else:
            windowStart = windowStarts[-1] + windowWidth
        windowStarts.append(windowStart)

    #convert the window starts and ends to megabases
    windowStarts = [x*1000000 for x in windowStarts]
    windowEnds = [x+(windowWidth*1000000)-1 for x in windowStarts]
    windowNames = methCNVDf.columns.values

    windowLocDf = pd.DataFrame(list(zip(windowNames,chromosomes,windowStarts,windowEnds)))
    windowLocDf.columns = ['windowName','chromosome','start','end']
    
    #get the location of all genes found in scRNA data (this is after merging
    #individual samples, so any patient gene position file suffices
    genePosFile = 'inferCNV_MGH122/MGH122_gene_pos_inferCNV.txt'
    genePosDf = pd.read_csv(genePosFile,index_col=0,header=None,sep='\t')
    genePosDf.columns = ['chromosome','start','end']
    genesInRNAFile = scRNACNVDf.columns.values
    genePosDf = genePosDf.loc[genesInRNAFile,:]
    
    #match each gene to a window
    geneLocList = list()
    for window in windowLocDf.iterrows():
        chromosome = window[1]['chromosome']
        start = window[1]['start']
        end = window[1]['end']
        genes = list(genePosDf.index[(genePosDf['chromosome'] == chromosome) & (((genePosDf['start'] >= start) & (genePosDf['start'] <= end)) | ((genePosDf['end'] >= start) & (genePosDf['end'] <= end)) | ((genePosDf['start'] <= start) & (genePosDf['end'] >= end)))])
        geneLocList.append(genes)
    windowLocDf['genes'] = geneLocList
    
    #get CNV for scRNA of all genes by window
    scRNACNVDfNew = pd.DataFrame()
    for window in windowLocDf.iterrows():
        scRNACNVDfNew[window[1].windowName] = scRNACNVDf.loc[:,window[1]['genes']].mean(axis=1)

    scRNACNVDfNew = scRNACNVDfNew.apply(np.log2)

    #clean the file with scRNA and scRRBs barcodes for matching
    mergingFile = 'scRNA_scRRBS_Barcodes_Matched.csv'
    mergingFileDf = pd.read_csv(mergingFile)
    mergingFileDf['scRNA_code'] = [mergingFileDf.iloc[x,0] + str(mergingFileDf.iloc[x,1]) for x in mergingFileDf.index]
    mergingFileDf['pool'] = [int(x[-1]) for x in mergingFileDf['Unnamed: 3']]
    mergingFileDf = mergingFileDf.drop(['scRNA','Unnamed: 1','Unnamed: 3'],axis=1)
    
    #extract the pool, barcode and sample from the index
    methCNVDf['scRRBS'] = [x.split('.')[1] for x in methCNVDf.index]
    methCNVDf['pool'] = [int(x.split('_')[-5][1]) for x in methCNVDf.index]
    methCNVDf['sample'] = [x.split('_')[1][3:] for x in methCNVDf.index]
    
    #extract scRNA_code and sample from the index
    scRNACNVDfNew.index = scRNACNVDfNew.index.str.replace('__', '_', regex=True)
    index = scRNACNVDfNew.index
    scRNACNVDfNew['sample'] = [x.split('_')[0][3:] for x in scRNACNVDfNew.index]
    scRNACNVDfNew['scRNA_code'] = [x.split('_')[-5] for x in scRNACNVDfNew.index]
    
    #merge the files
    scRNACNVDfNew = pd.merge(scRNACNVDfNew,mergingFileDf,how='left',left_on='scRNA_code',right_on='scRNA_code')
    
    scRNACNVDfNew.index = index
    scRNACNVDfNew = scRNACNVDfNew.reset_index().merge(methCNVDf,how='inner',
                                       left_on=['scRRBS','sample','pool'],
                                       right_on=['scRRBS','sample','pool'],
                                       suffixes=('_scRNA','_scRRBS')).set_index('index')
    
    #remove the now irrelevant keys
    scRNACNVDfNew = scRNACNVDfNew.drop(['sample','scRNA_code','scRRBS','pool'],axis=1)
    
    #add a column for the scRNA window without genes
    for i in range(len(scRNACNVDfNew.columns)):
        col = scRNACNVDfNew.columns[i]
        if (col[-7:] != '_scRRBS') & (col[-6:] != '_scRNA'):
            scRNACNVDfNew=scRNACNVDfNew.rename(columns = {col:col+'_scRRBS'})
            scRNACNVDfNew[col+'_scRNA'] = [np.nan]*scRNACNVDfNew.shape[0]
    
    columns = list(methCNVDf.columns.values)
    
    for x in ['scRRBS','pool','sample']:
        columns.remove(x)
    columnsx = [x+'_scRNA' for x in columns]
    columnsy = [x+'_scRRBS' for x in columns]
    columns = columnsx + columnsy
    scRNACNVDfNew = scRNACNVDfNew[columns]

    #save output for later
    scRNACNVDfNew.to_csv('CNV_w%s_s%s.csv' %(z,z))

#do the same for windows smaller than 1mb, note that these window sizes
#weren't used for further analysis
os.chdir('w0.1_s0.1')

for z in resamplesFiner:
    methCNVDf = methCNVFileImporter(samplesGBM,'_cpgs_w0.%s_s0.%s.csv' %(z,z))

    #filter out cells that did not pass QC
    methCNVFilterFile = 'allQC_DNAme_POST_filtered_50K_09112019.csv'
    methCNVFilterDf = pd.read_csv(methCNVFilterFile,index_col=0)
    methCNVFilter = list(methCNVFilterDf['Cell'])
    for i in range(len(methCNVFilter)):
        if methCNVFilter[i][:4] != 'GBM_':
            methCNVFilter[i] = 'GBM_' + methCNVFilter[i]
        
    for row in methCNVDf.iterrows():
        if '.'.join(row[0].split('.')[:2]) not in methCNVFilter:
            methCNVDf = methCNVDf.drop(row[0])
    
    methCNVDf.columns = [methCNVDf.columns.values[x].split('.')[0] + '.' + str(x+1) for x in range(len(methCNVDf.columns.values))]
        
    #import the scRNA CNV data
    samplesScRNA = ["105A","105B","105C",124,122]
    scRNADfList = list()
    #import the scRNA based CNV data
    for i in samplesScRNA:
        file = 'inferCNV_MGH%s/preliminary_modified_expression.csv' %i
        df = pd.read_csv(file,index_col=0)
        df = df.T
        annotationFile = 'inferCNV_MGH%s/MGH%s_annotation_file_inferCNV.txt' %(i,i)
        annoDf = pd.read_csv(annotationFile,sep='\t',index_col=0,header=None)
        df = df.join(annoDf)
        df = df[df.iloc[:,-1] == 'tumor']
        scRNADfList.append(df)
    
    #get the common genes in each sample
    common_cols = list(set.intersection(*(set(df.columns) for df in scRNADfList)))
    #concatenate using only the common columns
    scRNACNVDf = pd.concat([df[common_cols] for df in scRNADfList])
    scRNACNVDf = scRNACNVDf.drop(1,axis=1)
    
    #filter out cells that did not pass QC
    scRNAQCFilterFile = 'allQC_POST_filtered_scRNA_09122019.csv'
    scRNAQCFilterDf = pd.read_csv(scRNAQCFilterFile,index_col=0)
    scRNAQCFilter = list(scRNAQCFilterDf['X1'])
    for row in scRNACNVDf.iterrows():
        if row[0] not in scRNAQCFilter:
            scRNACNVDf = scRNACNVDf.drop(row[0])
    #find the window start and ends for all windows in the methylation file
    startsNew,endsNew = getChromosomeCutoffs(methCNVDf)
    windowWidth = z
    windowStarts = list()
    chromosomes = [int(x.split('.')[0][0:]) for x in methCNVDf.columns.values]
    
    for col in methCNVDf.columns.values:
        if col in startsNew:
            windowStart = 0
        else:
            windowStart = windowStarts[-1] + windowWidth
        windowStarts.append(windowStart)

    #convert the window starts and ends to megabases
    windowStarts = [x*1000000 for x in windowStarts]
    windowEnds = [x+(windowWidth*1000000)-1 for x in windowStarts]
    windowNames = methCNVDf.columns.values

    windowLocDf = pd.DataFrame(list(zip(windowNames,chromosomes,windowStarts,windowEnds)))
    windowLocDf.columns = ['windowName','chromosome','start','end']
    genePosFile = 'inferCNV_MGH122/MGH122_gene_pos_inferCNV.txt'
    genePosDf = pd.read_csv(genePosFile,index_col=0,header=None,sep='\t')
    genePosDf.columns = ['chromosome','start','end']
    genesInRNAFile = scRNACNVDf.columns.values
    genePosDf = genePosDf.loc[genesInRNAFile,:]
    
    geneLocList = list()
    for window in windowLocDf.iterrows():
        chromosome = window[1]['chromosome']
        start = window[1]['start']
        end = window[1]['end']
        genes = list(genePosDf.index[(genePosDf['chromosome'] == chromosome) & (((genePosDf['start'] >= start) & (genePosDf['start'] <= end)) | ((genePosDf['end'] >= start) & (genePosDf['end'] <= end)) | ((genePosDf['start'] <= start) & (genePosDf['end'] >= end)))])
        geneLocList.append(genes)
    windowLocDf['genes'] = geneLocList

    scRNACNVDfNew = pd.DataFrame()
    for window in windowLocDf.iterrows():
        scRNACNVDfNew[window[1].windowName] = scRNACNVDf.loc[:,window[1]['genes']].mean(axis=1)

    scRNACNVDfNew = scRNACNVDfNew.apply(np.log2)

    #get the merging file ready
    mergingFile = 'scRNA_scRRBS_Barcodes_Matched.csv'
    mergingFileDf = pd.read_csv(mergingFile)
    mergingFileDf['scRNA_code'] = [mergingFileDf.iloc[x,0] + str(mergingFileDf.iloc[x,1]) for x in mergingFileDf.index]
    mergingFileDf['pool'] = [int(x[-1]) for x in mergingFileDf['Unnamed: 3']]
    mergingFileDf = mergingFileDf.drop(['scRNA','Unnamed: 1','Unnamed: 3'],axis=1)
    
    #extract the pool, barcode and sample from the index
    methCNVDf['scRRBS'] = [x.split('.')[1] for x in methCNVDf.index]
    methCNVDf['pool'] = [int(x.split('_')[-5][1]) for x in methCNVDf.index]
    methCNVDf['sample'] = [x.split('_')[1][3:] for x in methCNVDf.index]
    
    #extract scRNA_code and sample from the index
    scRNACNVDfNew.index = scRNACNVDfNew.index.str.replace('__', '_', regex=True)
    index = scRNACNVDfNew.index
    scRNACNVDfNew['sample'] = [x.split('_')[0][3:] for x in scRNACNVDfNew.index]
    scRNACNVDfNew['scRNA_code'] = [x.split('_')[-5] for x in scRNACNVDfNew.index]
    
    #merge the merging file
    scRNACNVDfNew = pd.merge(scRNACNVDfNew,mergingFileDf,how='left',left_on='scRNA_code',right_on='scRNA_code')
    
    scRNACNVDfNew.index = index
    scRNACNVDfNew = scRNACNVDfNew.reset_index().merge(methCNVDf,how='inner',
                                       left_on=['scRRBS','sample','pool'],
                                       right_on=['scRRBS','sample','pool'],
                                       suffixes=('_scRNA','_scRRBS')).set_index('index')
    
    #remove the now irrelevant keys
    scRNACNVDfNew = scRNACNVDfNew.drop(['sample','scRNA_code','scRRBS','pool'],axis=1)
    
    #add a column for the scRNA window without genes
    for i in range(len(scRNACNVDfNew.columns)):
        col = scRNACNVDfNew.columns[i]
        if (col[-7:] != '_scRRBS') & (col[-6:] != '_scRNA'):
            scRNACNVDfNew=scRNACNVDfNew.rename(columns = {col:col+'_scRRBS'})
            scRNACNVDfNew[col+'_scRNA'] = [np.nan]*scRNACNVDfNew.shape[0]
    
    columns = list(methCNVDf.columns.values)
    
    for x in ['scRRBS','pool','sample']:
        columns.remove(x)
    columnsx = [x+'_scRNA' for x in columns]
    columnsy = [x+'_scRRBS' for x in columns]
    columns = columnsx + columnsy
    scRNACNVDfNew = scRNACNVDfNew[columns]

    scRNACNVDfNew.to_csv('CNV_w0.%s_s0.%s.csv' %(z,z))

#function to generate the quadrant plot with amplified/deleted predition
#allows for different thresholds to be classified as amplified or deleted
def visualizer(df,chrList,chr1,chr2,title,output,thresh):
    colNames = list()
    for i in chrList:
        curList = [x if (int(x.split('.')[0]) == i) & (x[-7:] == '_scRRBS') else '' for x in df.columns.values]
        curList = list(set(curList))
        curList.remove('')
        curList = [x.split('_')[0] for x in curList]
        colNames.extend(curList)
        
    scRNA = [str(x)+'_scRNA' for x in colNames]
    scRRBS = [str(x)+'_scRRBS' for x in colNames]

    chr1ColsScRRBS = [x if (int(x.split('.')[0]) == chr1) & (x[-7:] == '_scRRBS') else '' for x in df.columns.values]
    chr1ColsScRRBS = set(chr1ColsScRRBS)
    chr1ColsScRRBS.remove('')
    
    chr2ColsScRRBS = [x if (int(x.split('.')[0]) == chr2) & (x[-7:] == '_scRRBS') else '' for x in df.columns.values]
    chr2ColsScRRBS = set(chr2ColsScRRBS)
    chr2ColsScRRBS.remove('')

    chr1ColsScRNA = [x if (int(x.split('.')[0]) == chr1) & (x[-6:] == '_scRNA') else '' for x in df.columns.values]
    chr1ColsScRNA = set(chr1ColsScRNA)
    chr1ColsScRNA.remove('')
    
    chr2ColsScRNA = [x if (int(x.split('.')[0]) == chr2) & (x[-6:] == '_scRNA') else '' for x in df.columns.values]
    chr2ColsScRNA = set(chr2ColsScRNA)
    chr2ColsScRNA.remove('')
    
    plt.plot(figsize=[10,10])
    plt.axhline(linewidth=1, color='k', ls='--')
    plt.axvline(linewidth=1, color='k', ls='--')
    #for neutral, select all chromosomes that are not intersting and stack the df
    #to make the legend come out nice
    plt.scatter(x = df.loc[:,scRNA].stack(), y = df.loc[:,scRRBS].stack(),c='#D3D3D3',s=0.1, label='Neutral')
    plt.scatter(x = df.loc[:,list(chr1ColsScRNA)],y = df.loc[:,list(chr1ColsScRRBS)],c='r',s=0.1,label='Chr '+str(chr1))
    plt.scatter(x = df.loc[:,list(chr2ColsScRNA)],y = df.loc[:,list(chr2ColsScRRBS)],c='b',s=0.1, label='Chr '+str(chr2))
    plt.hlines([thresh,-thresh],[-1.5,-1.5],[1.5,1.5])
    plt.vlines([thresh,-thresh],[-1.5,-1.5],[1.5,1.5])
         
    plt.title(title)
    plt.xlabel('CNV based on scRNA')
    plt.ylabel('CNV based on scRRBS')
    plt.xlim([-1.5,1.5])
    plt.ylim([-1.5,1.5])
    plt.yticks([-1,0,1])
    plt.xticks([-1,0,1])
    lgnd = plt.legend()
    #change the marker size manually for both lines
    lgnd.legendHandles[0]._sizes = [30]
    lgnd.legendHandles[1]._sizes = [30]
    lgnd.legendHandles[2]._sizes = [30]
    
    plt.savefig(output)
    plt.show()

normalChromosomes = [1,2]
amplifiedChromosomes = [7]
deletedChromosomes = [10]

#this function was run using window size = 20mb
visualizer(scRNACNVDfNew,normalChromosomes,amplifiedChromosomes[0],deletedChromosomes[0],'Comparison of scRRBS and scRNA inference of CNV \n Threshold = 0.15','Comparison_of_scRRBS_and_scRNA_inference_of_CNV_chromosome_sections_0.15.pdf',0.15)
visualizer(scRNACNVDfNew,normalChromosomes,amplifiedChromosomes[0],deletedChromosomes[0],'Comparison of scRRBS and scRNA inference of CNV \n Threshold = 0.3','Comparison_of_scRRBS_and_scRNA_inference_of_CNV_chromosome_sections_0.3.pdf',0.30)

#get the windows in the chromosomes of interest
def colForROC(df,chromosomes):
    currentDf = df.iloc[:,int(df.shape[1]/2):]
    currentDf = currentDf.dropna(axis=1)
    colsList = list()

    for col in currentDf.columns.values:
        if int(col.split('.')[0]) in chromosomes:
           colsList.append(col)
           colsList.append(col[:-7]+'_scRNA')
    return(colsList)
    
#for each window, say if it should be amplified, deleted or neutral based on chromosome
def conditionsForROC(chromosomes,conditions,columns):
    conditionsOutputList = list()
    for col in columns:
        condition = conditions[chromosomes.index(int(col.split('.')[0]))]
        conditionsOutputList.append(condition)
    return(conditionsOutputList)

#split the dataframe into the scRNA and scRRBS
#scRNA data needs to be weighted by number of genes per segment
def chrExtractor(df,chromosome):
    columns = colForROC(df,chromosome)
    columnsScRNA = list()
    columnsscRRBS = list()
    for col in columns:
        if col[-7:] == '_scRRBS':
            columnsscRRBS.append(col)
        if col[-6:] == '_scRNA':
            columnsScRNA.append(col)
    return([columnsscRRBS,columnsScRNA])

#finds the confidence interval cutoffs for the ROCPlotter
def cutoffFinder(inputList):
    lowerValues = list()
    upperValues = list()
    for i in inputList:
        i.sort()
        lowerValue = i[int(round(0.05 * len(i),0))]
        upperValue = i[int(round(0.95 * len(i),0))]
        lowerValues.append(lowerValue)
        upperValues.append(upperValue)
    return(lowerValues,upperValues)
    
#finds the confidence interval cutoffs for AUC used in the ROCPlotter
def cutoffFinderAuc(inputList):
    inputList.sort()
    lowerValue = inputList[int(round(0.05 * len(inputList),0))]
    upperValue = inputList[int(round(0.95 * len(inputList),0))]

    return(lowerValue,upperValue)
 
#plots ROC curves with confidence intervals generated by bootstrapping
def ROCPlotter(df,cols,conditions,thresholdList,posNeg,title,output):
    columns = colForROC(df,cols)
    columnsScRNA = list()
    columnsScRRBS = list()
    
    for col in columns:
        if col[-6:] == '_scRNA':
            columnsScRNA.append(col)
        elif col[-7:] == '_scRRBS':
            columnsScRRBS.append(col)
            
    conditionsScRNA = conditionsForROC(cols,conditions,columnsScRNA)
    conditionsScRRBS = conditionsForROC(cols,conditions,columnsScRRBS)
    
    tprscRNAList = []
    fprscRNAList = []
    tprscRRBSList = []
    fprscRRBSList = []
    
    n_bootstraps = 200
    
    scRNATPRList = list()
    scRNAFPRList = list()
    
    scRRBSTPRList = list()
    scRRBSFPRList = list()
    
    for thresh in thresholdList:
        print('threshold = %s' %thresh)
        scRNAyTrue, scRNAyPred = binaryROCInputGenerator(df,columnsScRNA,conditionsScRNA,thresh,posNeg)
        matrix = confusion_matrix(scRNAyTrue, scRNAyPred)
        try:
            tn, fp, fn, tp = matrix.ravel()
            tpr = tp/(tp+fn)
            fpr = fp/(fp+tn)
            tprscRNAList.append(tpr)
            fprscRNAList.append(fpr)
        except:
            print('insufficient data')
        
        
        scRRBSyTrue, scRRBSyPred = binaryROCInputGenerator(df,columnsScRRBS,conditionsScRRBS,thresh,posNeg)
        matrix = confusion_matrix(scRRBSyTrue, scRRBSyPred)
        try:
            tn, fp, fn, tp = matrix.ravel()
            tpr = tp/(tp+fn)
            fpr = fp/(fp+tn)
            tprscRRBSList.append(tpr)
            fprscRRBSList.append(fpr)
        except:
            print('insufficient data')
            
        scRNATPR = list()
        scRNAFPR = list()
        
        scRRBSTPR = list()
        scRRBSFPR = list()
        
        #bootstraps should be run on whole cells, so the data transformation has to be rerun
        for i in range(n_bootstraps):
            bootstrapDf = resample(df,n_samples=int(round(scRNACNVDfNew.shape[0]*0.6,0)))
    
            bootstrapColumns = colForROC(bootstrapDf,cols)
            bootstrapColumnsScRNA = list()
            bootstrapColumnsScRRBS = list()
            
            for col in bootstrapColumns:
                if col[-6:] == '_scRNA':
                    bootstrapColumnsScRNA.append(col)
                elif col[-7:] == '_scRRBS':
                    bootstrapColumnsScRRBS.append(col)
            bootstrapConditionsScRNA = conditionsForROC(cols,conditions,bootstrapColumnsScRNA)
            bootstrapConditionsScRRBS = conditionsForROC(cols,conditions,bootstrapColumnsScRRBS)   
    
            bootstrapScRNAyTrue, bootstrapScRNAyPred = binaryROCInputGenerator(bootstrapDf,bootstrapColumnsScRNA,bootstrapConditionsScRNA,thresh,posNeg)
            bootstrapScRRBSyTrue, bootstrapScRRBSyPred = binaryROCInputGenerator(bootstrapDf,bootstrapColumnsScRRBS,bootstrapConditionsScRRBS,thresh,posNeg)
            
            matrix = confusion_matrix(bootstrapScRNAyTrue, bootstrapScRNAyPred)  
            try:
                tn, fp, fn, tp = matrix.ravel()
                tpr = tp/(tp+fn)
                fpr = fp/(fp+tn)
                scRNATPR.append(tpr)
                scRNAFPR.append(fpr)
            except:
                pass
                
            matrix = confusion_matrix(bootstrapScRRBSyTrue, bootstrapScRRBSyPred)
            try:
                tn, fp, fn, tp = matrix.ravel()
                tpr = tp/(tp+fn)
                fpr = fp/(fp+tn)
                scRRBSTPR.append(tpr)
                scRRBSFPR.append(fpr)
            except:
                pass
                
        scRNATPRList.append(scRNATPR)
        scRNAFPRList.append(scRNAFPR)
    
        scRRBSTPRList.append(scRRBSTPR)
        scRRBSFPRList.append(scRRBSFPR)           
    
    scRRBSAuc = auc(fprscRRBSList,tprscRRBSList)
    scRNAAuc = auc(fprscRNAList,tprscRNAList)
    
    scRNATPRLow,scRNATPRUp = cutoffFinder(scRNATPRList)
    
    scRRBSTPRLow,scRRBSTPRUp = cutoffFinder(scRRBSTPRList)
     
    scRRBSAuc = auc(fprscRRBSList,tprscRRBSList)
    scRNAAuc = auc(fprscRNAList,tprscRNAList)
    plt.plot()
    plt.plot(fprscRNAList,tprscRNAList,label='scRNA AUC=%s' %round(scRNAAuc,2))
    plt.plot(fprscRRBSList,tprscRRBSList,label='scRRBS AUC=%s' %round(scRRBSAuc,2))
    plt.fill_between(fprscRNAList,scRNATPRLow,scRNATPRUp,alpha=0.5)
    plt.fill_between(fprscRRBSList,scRRBSTPRLow,scRRBSTPRUp,alpha=0.5)
    plt.title(title)
    plt.xlabel('False negative rate')
    plt.ylabel('True positive rate')
    plt.legend()
    plt.savefig(output)
    plt.show()

#converts amplification/deletion/neutral into a binary classification result
def binaryROCInputGenerator(df,columns,conditions,thresh,posNeg):
    yTrue = list()
    for i in range(len(conditions)):
        condition = conditions[i]
        if posNeg == 'pos':
            if condition == ">":
                interimYTrue = [1]*df.shape[0]
                yTrue.extend(interimYTrue)
            else:
                interimYTrue = [-1]*df.shape[0]
                yTrue.extend(interimYTrue)
        else:
            if condition == "<":
                interimYTrue = [1]*df.shape[0]
                yTrue.extend(interimYTrue)
            else:
                interimYTrue = [-1]*df.shape[0]
                yTrue.extend(interimYTrue)            
    yPred = list()
    for i in range(len(columns)):
        column = columns[i]
        if posNeg == 'pos':
            interimYPred = [1 if x > thresh else -1 for x in df[column]]
        else:
            interimYPred = [1 if x < thresh else -1 for x in df[column]]
        yPred.extend(interimYPred)
    return(yTrue,yPred)

#gets the AUC scores with confidence interval generated through bootstrapping
def AUCScores(df,cols,conditions,thresholdList,posNeg):
    columns = colForROC(df,cols)
    columnsScRNA = list()
    columnsScRRBS = list()

    for col in columns:
        if col[-6:] == '_scRNA':
            columnsScRNA.append(col)
        elif col[-7:] == '_scRRBS':
            columnsScRRBS.append(col)
           
    conditionsScRNA = conditionsForROC(cols,conditions,columnsScRNA)
    conditionsScRRBS = conditionsForROC(cols,conditions,columnsScRRBS)
    
    tprscRNAList = []
    fprscRNAList = []
    tprscRRBSList = []
    fprscRRBSList = []
    
    n_bootstraps = 200

    for thresh in thresholdList:
        print('thresh = %s' %thresh)
        scRNAyTrue, scRNAyPred = binaryROCInputGenerator(df,columnsScRNA,conditionsScRNA,thresh,posNeg)
        matrix = confusion_matrix(scRNAyTrue, scRNAyPred)
        try:
            tn, fp, fn, tp = matrix.ravel()
            tpr = tp/(tp+fn)
            fpr = fp/(fp+tn)
            tprscRNAList.append(tpr)
            fprscRNAList.append(fpr)
        except:
            print('insufficient data')
        
        
        scRRBSyTrue, scRRBSyPred = binaryROCInputGenerator(df,columnsScRRBS,conditionsScRRBS,thresh,posNeg)
        matrix = confusion_matrix(scRRBSyTrue, scRRBSyPred)  
        try:
            tn, fp, fn, tp = matrix.ravel()
            tpr = tp/(tp+fn)
            fpr = fp/(fp+tn)
            tprscRRBSList.append(tpr)
            fprscRRBSList.append(fpr)
        except:
            print('insufficient data')

    scRNAAucBootList = list()
    scRRBSAucBootList = list()
    
    scRNATPRBootList = list()
    scRNAFPRBootList = list()
    scRRBSTPRBootList = list()
    scRRBSFPRBootList = list()
    
    #bootstraps should be run on whole cells, so the data transformation has to be rerun
    for i in range(n_bootstraps):
        print('bootstrap = %s' %i)
        bootstrapDf = resample(df,n_samples=int(round(df.shape[0]*0.6,0)),random_state = i)

        bootstrapColumns = colForROC(bootstrapDf,cols)
        bootstrapColumnsScRNA = list()
        bootstrapColumnsScRRBS = list()
        
        for col in bootstrapColumns:
            if col[-6:] == '_scRNA':
                bootstrapColumnsScRNA.append(col)
            elif col[-7:] == '_scRRBS':
                bootstrapColumnsScRRBS.append(col)
        
        bootstrapConditionsScRNA = conditionsForROC(cols,conditions,bootstrapColumnsScRNA)
        bootstrapConditionsScRRBS = conditionsForROC(cols,conditions,bootstrapColumnsScRRBS)   
        
            
        scRNATPR = list()
        scRNAFPR = list()
        
        scRRBSTPR = list()
        scRRBSFPR = list()
        
        for thresh in thresholdList:
            bootstrapScRNAyTrue, bootstrapScRNAyPred = binaryROCInputGenerator(bootstrapDf,bootstrapColumnsScRNA,bootstrapConditionsScRNA,thresh,posNeg)
            bootstrapScRRBSyTrue, bootstrapScRRBSyPred = binaryROCInputGenerator(bootstrapDf,bootstrapColumnsScRRBS,bootstrapConditionsScRRBS,thresh,posNeg)
        
            matrix = confusion_matrix(bootstrapScRNAyTrue, bootstrapScRNAyPred)  
            try:
                tn, fp, fn, tp = matrix.ravel()
                tpr = tp/(tp+fn)
                fpr = fp/(fp+tn)
                scRNATPR.append(tpr)
                scRNAFPR.append(fpr)
            except:
                print('insufficient data')
                
            matrix = confusion_matrix(bootstrapScRRBSyTrue, bootstrapScRRBSyPred)  
            try:
                tn, fp, fn, tp = matrix.ravel()
                tpr = tp/(tp+fn)
                fpr = fp/(fp+tn)
                scRRBSTPR.append(tpr)
                scRRBSFPR.append(fpr)
            except:
                print('insufficient data')
        
        scRNAAucBoot = auc(scRNAFPR,scRNATPR)
        scRNAAucBootList.append(scRNAAucBoot)

        scRRBSAucBoot = auc(scRRBSFPR,scRRBSTPR)
        scRRBSAucBootList.append(scRRBSAucBoot)
        
        scRNATPRBootList.append(scRNATPR)
        scRNAFPRBootList.append(scRNAFPR)
        scRRBSTPRBootList.append(scRRBSTPR)
        scRRBSFPRBootList.append(scRRBSFPR)
        
    scRNAAuc = auc(fprscRNAList,tprscRNAList)
   
    scRRBSAuc = auc(fprscRRBSList,tprscRRBSList)


    scRNAAucLow,scRNAAucUp = cutoffFinderAuc(scRNAAucBootList)
    
    scRRBSAucLow,scRRBSAucUp = cutoffFinderAuc(scRRBSAucBootList)
    
    outputList = [scRNAAuc,scRRBSAuc,scRNAAucLow,scRNAAucUp,scRRBSAucLow,scRRBSAucUp]
    rawDataList = [fprscRNAList,tprscRNAList,fprscRRBSList,tprscRRBSList]

    return(outputList,rawDataList,scRNATPRBootList,scRNAFPRBootList,scRRBSTPRBootList,scRRBSFPRBootList)

#downsamples scRNA data for different window sizes
def downSample(df,resample):
    cols = df.columns.values[::resample]
    for i in range(len(cols)):
        if i < len(cols)-1:
            df[cols[i]] = df.loc[:,cols[i]:df.columns.values[list(df.columns.values).index(cols[i])+resample-1]].mean(axis=1)
        else:
            df[cols[i]] = df.loc[:,cols[i]:].mean(axis=1)
    df = df[cols]
    return(df)
    
#downsamples scRNA data for different window sizes
def downSampleScRNA(df,resample,geneDf):
    numGenesPerWindow = geneDf.to_frame()
    numGenesPerWindow.columns = ['chromosome']
    for i in df.columns.values:
        try:
            df.loc[:,i] = df.loc[:,i] * numGenesPerWindow.loc[i[:-6],'chromosome']
        except:
            df.loc[:,i] = df.loc[:,i]
    cols = df.columns.values[::resample]
    for i in range(len(cols)):
        if i < len(cols)-1:
            numGenes = numGenesPerWindow['chromosome'].to_frame().loc[cols[i][:-6]:df.columns.values[list(df.columns.values).index(cols[i])+resample-1][:-6]].sum(axis=0)
            df[cols[i]] = df.loc[:,cols[i]:df.columns.values[list(df.columns.values).index(cols[i])+resample-1]].sum(axis=1)/numGenes[0]
        else:
            numGenes = numGenesPerWindow['chromosome'].to_frame().loc[cols[i][:-6]:].sum(axis=0)
            df[cols[i]] = df.loc[:,cols[i]:].mean(axis=1)/numGenes[0]
    df = df[cols]
    return(df)

#get all predictions, etc. for different thresholds
resamples = [1,2,3,5,10,20,30,40,50]
thresholds = [x/100 for x in range(-100,101)]
columns = [1,2,7,10]
conditions = ['neu','neu','>','<']

aucDict = dict()
for z in resamples:
    start_time = time.time()
    chr1scRNA = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[1])[1]]
    chr2scRNA = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[2])[1]]
    chr7scRNA = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[7])[1]]
    chr10scRNA = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[10])[1]]

    chr1scRNAResampled = downSampleScRNA(chr1scRNA,z,genesPerLocDf)
    chr2scRNAResampled = downSampleScRNA(chr2scRNA,z,genesPerLocDf)
    chr7scRNAResampled = downSampleScRNA(chr7scRNA,z,genesPerLocDf)
    chr10scRNAResampled = downSampleScRNA(chr10scRNA,z,genesPerLocDf) 

    chr1scRRBS = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[1])[0]]
    chr2scRRBS = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[2])[0]]
    chr7scRRBS = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[7])[0]]
    chr10scRRBS = scRNACNVDfNew[chrExtractor(scRNACNVDfNew,[10])[0]]

    chr1scRRBSResampled = downSample(chr1scRRBS,z)
    chr2scRRBSResampled = downSample(chr2scRRBS,z)
    chr7scRRBSResampled = downSample(chr7scRRBS,z)
    chr10scRRBSResampled = downSample(chr10scRRBS,z) 

    dfList = [chr1scRNAResampled,chr2scRNAResampled,chr7scRNAResampled,chr10scRNAResampled,chr1scRRBSResampled,chr2scRRBSResampled,chr7scRRBSResampled,chr10scRRBSResampled]
    
    #merge the dataframes back into one
    workingDf = pd.concat(dfList,axis=1)

    output,rawDataList,scRNATPRBootList,scRNAFPRBootList,scRRBSTPRBootList,scRRBSFPRBootList = AUCScores(workingDf,columns,conditions,thresholds,'pos')
    
    aucDict['%smb' %z] = output
    
    #all these intermediary outputs are saved for posterity but not used
    rawDataDf = pd.DataFrame(rawDataList)
    rawDataDf = rawDataDf.T
    rawDataDf.columns = ['scRNA fpr','scRNA tpr','scRRBS fpr','scRRBS tpr']
    rawDataDf.to_csv('%smb_non_bootstrap_output.csv' %z)
    
    scRNATPRBootDf = pd.DataFrame(scRNATPRBootList)
    scRNATPRBootDf = scRNATPRBootDf.T
    scRNATPRBootDf.columns = ['bootstrap '+str(x) for x in range(200)]
    scRNATPRBootDf.to_csv('%smb_bootstrap_scRNA_tpr_output.csv' %z)
    
    scRNAFPRBootDf = pd.DataFrame(scRNAFPRBootList)
    scRNAFPRBootDf = scRNAFPRBootDf.T
    scRNAFPRBootDf.columns = ['bootstrap '+str(x) for x in range(200)]
    scRNAFPRBootDf.to_csv('%smb_bootstrap_scRNA_fpr_output.csv' %z)
    
    scRRBSTPRBootDf = pd.DataFrame(scRRBSTPRBootList)
    scRRBSTPRBootDf = scRRBSTPRBootDf.T
    scRRBSTPRBootDf.columns = ['bootstrap '+str(x) for x in range(200)]
    scRRBSTPRBootDf.to_csv('%smb_bootstrap_scRRBS_tpr_output.csv' %z)
    
    scRRNSFPRBootDf = pd.DataFrame(scRRBSFPRBootList)
    scRRNSFPRBootDf = scRRNSFPRBootDf.T
    scRRNSFPRBootDf.columns = ['bootstrap '+str(x) for x in range(200)]
    scRRNSFPRBootDf.to_csv('%smb_bootstrap_scRRBS_fpr_output.csv' %z)
    print("--- %s seconds for resample %s---" % ((time.time() - start_time),z))

with open('aucDict.csv', 'w') as f:
    w = csv.DictWriter(f, aucDict.keys())
    w.writeheader()
    w.writerow(aucDict)
    
#function to better navigate the AUCdict
def dictUnpacker(dictionary,partialKeys,item):
    output = [dictionary['%smb' %i][item] for i in partialKeys]
    return output

resamples = [0.1,0.5,1,2,3,5,10,20,30,40,50]
scRNAAucForPlot = dictUnpacker(aucDict,resamples,0)
scRRBSAucForPlot = dictUnpacker(aucDict,resamples,1)
scRNAAucLowForPlot = dictUnpacker(aucDict,resamples,2)
scRNAAucUpForPlot = dictUnpacker(aucDict,resamples,3)
scRRBSAucLowForPlot = dictUnpacker(aucDict,resamples,4)
scRRBSAucUpForPlot = dictUnpacker(aucDict,resamples,5)

#this is used to generate Fig 1d bottom (to generate supp. Fig 3c, change 'pos' to 'neg')
resamplesROC = [1,5,10,20,30,40,50]
thresholds = [x/100 for x in range(-100,101)]
columns = [1,2,7,10]
conditions = ['neu','neu','>','<']
for z in resamplesROC:
    start_time = time.time()
    
    workingDf = pd.read_csv('CNV_w%s_s%s.csv' %(z,z),index_col=0)
    ROCPlotter(workingDf,columns,conditions,thresholds,'pos','ROC of CNVs based on scRNA/scRRBS \n Amplified chromosomes vs non-amplified %smb' %z,'ROC_Amplified_vs_non_amplified_chromosome_%smb_sections.pdf' %z)
    print("--- %s seconds for resample %s---" % ((time.time() - start_time),z))

#generate the files for whole chromosome ROC. Note: whole chromosome
#analysis was done separately, so the code is a copy of the above
os.chdir('w1_s1')

methCNVDf = methCNVFileImporter(samplesGBM,'_cpgs_whole_chromosomes.csv')

#filter out cells that did not pass QC
methCNVFilterFile = 'allQC_DNAme_POST_filtered_50K_09112019.csv'
methCNVFilterDf = pd.read_csv(methCNVFilterFile,index_col=0)
methCNVFilter = list(methCNVFilterDf['Cell'])
for i in range(len(methCNVFilter)):
    if methCNVFilter[i][:4] != 'GBM_':
        methCNVFilter[i] = 'GBM_' + methCNVFilter[i]
    
for row in methCNVDf.iterrows():
    if '.'.join(row[0].split('.')[:2]) not in methCNVFilter:
        methCNVDf = methCNVDf.drop(row[0])

methCNVDf.columns = [methCNVDf.columns.values[x].split('.')[0] + '.' + str(x+1) for x in range(len(methCNVDf.columns.values))]
    
#import the scRNA CNV data
samplesScRNA = ["105A","105B","105C",124,122]
scRNADfList = list()
#import the scRNA based CNV data
for i in samplesScRNA:
    file = 'inferCNV_MGH%s/preliminary_modified_expression.csv' %i
    df = pd.read_csv(file,index_col=0)
    df = df.T
    annotationFile = 'inferCNV_MGH%s/MGH%s_annotation_file_inferCNV.txt' %(i,i)
    annoDf = pd.read_csv(annotationFile,sep='\t',index_col=0,header=None)
    df = df.join(annoDf)
    df = df[df.iloc[:,-1] == 'tumor']
    scRNADfList.append(df)

#get the common genes in each sample
common_cols = list(set.intersection(*(set(df.columns) for df in scRNADfList)))
#concatenate using only the common columns
scRNACNVDf = pd.concat([df[common_cols] for df in scRNADfList])
scRNACNVDf = scRNACNVDf.drop(1,axis=1)

#filter out cells that did not pass QC
scRNAQCFilterFile = 'allQC_POST_filtered_scRNA_09122019.csv'
scRNAQCFilterDf = pd.read_csv(scRNAQCFilterFile,index_col=0)
scRNAQCFilter = list(scRNAQCFilterDf['X1'])
for row in scRNACNVDf.iterrows():
    if row[0] not in scRNAQCFilter:
        scRNACNVDf = scRNACNVDf.drop(row[0])
#find the window start and ends for all windows in the methylation file
startsNew,endsNew = getChromosomeCutoffs(methCNVDf)
windowWidth = z
windowStarts = list()
chromosomes = [int(x.split('.')[0][0:]) for x in methCNVDf.columns.values]

for col in methCNVDf.columns.values:
    if col in startsNew:
        windowStart = 0
    else:
        windowStart = windowStarts[-1] + windowWidth
    windowStarts.append(windowStart)

#convert the window starts and ends to megabases
windowStarts = [x*1000000 for x in windowStarts]
windowEnds = [x+(windowWidth*1000000)-1 for x in windowStarts]
windowNames = methCNVDf.columns.values

windowLocDf = pd.DataFrame(list(zip(windowNames,chromosomes,windowStarts,windowEnds)))
windowLocDf.columns = ['windowName','chromosome','start','end']
genePosFile = 'inferCNV_MGH122/MGH122_gene_pos_inferCNV.txt'
genePosDf = pd.read_csv(genePosFile,index_col=0,header=None,sep='\t')
genePosDf.columns = ['chromosome','start','end']
genesInRNAFile = scRNACNVDf.columns.values
genePosDf = genePosDf.loc[genesInRNAFile,:]

geneLocList = list()
for window in windowLocDf.iterrows():
    chromosome = window[1]['chromosome']
    start = window[1]['start']
    end = window[1]['end']
    genes = list(genePosDf.index[(genePosDf['chromosome'] == chromosome) & (((genePosDf['start'] >= start) & (genePosDf['start'] <= end)) | ((genePosDf['end'] >= start) & (genePosDf['end'] <= end)) | ((genePosDf['start'] <= start) & (genePosDf['end'] >= end)))])
    geneLocList.append(genes)
windowLocDf['genes'] = geneLocList

scRNACNVDfNew = pd.DataFrame()
for window in windowLocDf.iterrows():
    scRNACNVDfNew[window[1].windowName] = scRNACNVDf.loc[:,window[1]['genes']].mean(axis=1)

scRNACNVDfNew = scRNACNVDfNew.apply(np.log2)

#get the merging file ready
mergingFile = 'GBM/scRNA_scRRBS_Barcodes_Matched.csv'
mergingFileDf = pd.read_csv(mergingFile)
mergingFileDf['scRNA_code'] = [mergingFileDf.iloc[x,0] + str(mergingFileDf.iloc[x,1]) for x in mergingFileDf.index]
mergingFileDf['pool'] = [int(x[-1]) for x in mergingFileDf['Unnamed: 3']]
mergingFileDf = mergingFileDf.drop(['scRNA','Unnamed: 1','Unnamed: 3'],axis=1)

#extract the pool, barcode and sample from the index
methCNVDf['scRRBS'] = [x.split('.')[1] for x in methCNVDf.index]
methCNVDf['pool'] = [int(x.split('_')[-5][1]) for x in methCNVDf.index]
methCNVDf['sample'] = [x.split('_')[1][3:] for x in methCNVDf.index]

#extract scRNA_code and sample from the index
scRNACNVDfNew.index = scRNACNVDfNew.index.str.replace('__', '_', regex=True)
index = scRNACNVDfNew.index
scRNACNVDfNew['sample'] = [x.split('_')[0][3:] for x in scRNACNVDfNew.index]
scRNACNVDfNew['scRNA_code'] = [x.split('_')[-5] for x in scRNACNVDfNew.index]

#merge the merging file
scRNACNVDfNew = pd.merge(scRNACNVDfNew,mergingFileDf,how='left',left_on='scRNA_code',right_on='scRNA_code')

scRNACNVDfNew.index = index
scRNACNVDfNew = scRNACNVDfNew.reset_index().merge(methCNVDf,how='inner',
                                   left_on=['scRRBS','sample','pool'],
                                   right_on=['scRRBS','sample','pool'],
                                   suffixes=('_scRNA','_scRRBS')).set_index('index')

#remove the now irrelevant keys
scRNACNVDfNew = scRNACNVDfNew.drop(['sample','scRNA_code','scRRBS','pool'],axis=1)

#add a column for the scRNA window without genes
for i in range(len(scRNACNVDfNew.columns)):
    col = scRNACNVDfNew.columns[i]
    if (col[-7:] != '_scRRBS') & (col[-6:] != '_scRNA'):
        scRNACNVDfNew=scRNACNVDfNew.rename(columns = {col:col+'_scRRBS'})
        scRNACNVDfNew[col+'_scRNA'] = [np.nan]*scRNACNVDfNew.shape[0]

columns = list(methCNVDf.columns.values)

for x in ['scRRBS','pool','sample']:
    columns.remove(x)
columnsx = [x+'_scRNA' for x in columns]
columnsy = [x+'_scRRBS' for x in columns]
columns = columnsx + columnsy
scRNACNVDfNew = scRNACNVDfNew[columns]

scRNACNVDfNew.to_csv('CNV_whole_chromosomes.csv')

#generate the whole chromosome ROC
thresholds = [x/100 for x in range(-100,101)]
columns = [1,2,7,10]
conditions = ['neu','neu','>','<']
workingDf = pd.read_csv('CNV_whole_chromosomes.csv',index_col=0)
ROCPlotter(workingDf,columns,conditions,thresholds,'pos','ROC of CNVs based on scRNA/scRRBS \n Amplified chromosomes vs non-amplified whole chromosomes','ROC_Amplified_vs_non_amplified_chromosome_whole_chromosomes.pdf')

aucDict = dict()

output,rawDataList,scRNATPRBootList,scRNAFPRBootList,scRRBSTPRBootList,scRRBSFPRBootList = AUCScores(workingDf,columns,conditions,thresholds,'pos')

aucDict['whole_chromosome'] = output
with open('/Users/lkluegel/Desktop/auc_output2/aucDict_amplified_whole_chromosomes.csv', 'w') as f:  
    w = csv.DictWriter(f, aucDict.keys())
    w.writeheader()
    w.writerow(aucDict)
    
#AUC data with whole chromosome data
reader = csv.DictReader(open("aucDict_amplified.csv"))
dictobj = next(reader) 
aucDict = dict()
for key in dictobj.keys():
    aucDict[key] = dictobj[key]

#remove small windows as they aren't valuable
del aucDict['0.1mb']
del aucDict['0.5mb']

for k,v in aucDict.items():
    item = v[1:-1].split(',')
    item = [float(x) for x in item]
    aucDict[k] = item


    
#generate AUC plot (Fig 1d top) for windows of interest, pulling out data of interest
resamples = [1,2,3,5,10,20,30,40,50]
scRNAAucForPlot = dictUnpacker(aucDict,resamples,0)
scRRBSAucForPlot = dictUnpacker(aucDict,resamples,1)
scRNAAucLowForPlot = dictUnpacker(aucDict,resamples,2)
scRNAAucUpForPlot = dictUnpacker(aucDict,resamples,3)
scRRBSAucLowForPlot = dictUnpacker(aucDict,resamples,4)
scRRBSAucUpForPlot = dictUnpacker(aucDict,resamples,5)

#get the whole chromosome data
reader = csv.DictReader(open("aucDict_amplified_whole_chromosomes.csv"))
dictobj = next(reader) 
wholeChr = dictobj['whole_chromosome']
wholeChr = wholeChr[1:-1].split(',')
wholeChr = [float(x) for x in wholeChr]

scRNAAucForPlot.append(wholeChr[0])
scRRBSAucForPlot.append(wholeChr[1])
scRNAAucLowForPlot.append(wholeChr[2])
scRNAAucUpForPlot.append(wholeChr[3])
scRRBSAucLowForPlot.append(wholeChr[4])
scRRBSAucUpForPlot.append(wholeChr[5])

#function to smooth AUC curve
def fivepl(x, a, b, c, d, g):
    return(((a-d)/((1+((x/c)** b))**g))+d)
    
#smooth the data
popt1, pcov1 = curve_fit(fivepl, resamples, scRNAAucForPlot[:-1])
popt2, pcov2 = curve_fit(fivepl, resamples, scRRBSAucForPlot[:-1])
popt3, pcov3 = curve_fit(fivepl, resamples, scRNAAucLowForPlot[:-1])
popt4, pcov4 = curve_fit(fivepl, resamples, scRNAAucUpForPlot[:-1])
popt5, pcov5 = curve_fit(fivepl, resamples, scRRBSAucLowForPlot[:-1])
popt6, pcov6 = curve_fit(fivepl, resamples, scRRBSAucUpForPlot[:-1])

xs = np.linspace(1,50,100)

plt.plot()
plt.plot(xs,fivepl(xs,*popt1),label='scRNA')
plt.plot(xs,fivepl(xs,*popt2),label='scDNAme')
plt.plot([50,70],[fivepl(50,*popt1),scRNAAucForPlot[-1]],linestyle='dashed',c='#1f77b4')
plt.plot([50,70],[fivepl(50,*popt2),scRRBSAucForPlot[-1]],linestyle='dashed',c='#ff7f0e')
plt.fill_between(xs,fivepl(xs,*popt3),fivepl(xs,*popt4),alpha=0.5)
plt.fill_between(xs,fivepl(xs,*popt5),fivepl(xs,*popt6),alpha=0.5)
plt.legend()
plt.xlabel('Chromosome segment size in mb')
plt.ylabel('AUC for prediction of \n amplification vs non-amplification')
plt.title('Smoothed AUC for prediction of amplification vs non-amplification \n by chromosome segment size scDNAme vs scRNA')
plt.xticks((1,5,10,20,30,40,50,70),[1,5,10,20,30,40,50,'whole \n chr'])
plt.savefig('AUC_smoothed_amplified_with_whole_chromosomes.pdf')
plt.show()

#find correlation between scRRBS and scRNA based CNV for Fig 1c
file = 'CNV_w20_s20.csv'
workingDf = pd.read_csv(file,index_col=0)
columns = colForROC(workingDf,[1,2,7,10])
columnsScRNA = list()
columnsScRRBS = list()

for col in columns:
    if col[-6:] == '_scRNA':
        columnsScRNA.append(col)
    elif col[-7:] == '_scRRBS':
        columnsScRRBS.append(col)

X = workingDf.loc[:,columnsScRNA].stack()
Y = workingDf.loc[:,columnsScRRBS].stack()

correlation = pearsonr(X,Y)
correlation = correlation[0]


def r_to_z(r):
    return math.log((1 + r) / (1 - r)) / 2.0

def z_to_r(z):
    e = math.exp(2 * z)
    return((e - 1) / (e + 1))

def r_confidence_interval(r, alpha, n):
    z = r_to_z(r)
    se = 1.0 / math.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha/2)  # 2-tailed z critical value

    lo = z - z_crit * se
    hi = z + z_crit * se

    # Return a sequence
    return (z_to_r(lo), z_to_r(hi))

z = r_to_z(correlation)
normalizedCorrelation = z_to_r(z)
low,high = r_confidence_interval(normalizedCorrelation,0.95,len(X))
print(low)
print(correlation)
print(high)
