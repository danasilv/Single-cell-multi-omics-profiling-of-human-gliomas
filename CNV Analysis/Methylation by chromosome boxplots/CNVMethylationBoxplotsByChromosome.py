import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations
import numpy as np
import os

#finds chromosome cutoffs in datafiles
def getChromosomeCutoffs(df):
    starts = list()
    ends = list()
    columnNames = df.columns
    a = [x.split(".") for x in columnNames]
    
    for i in range(len(a)-1):
        rightChr = int(a[i+1][0][3:])
        leftChr = int(a[i][0][3:])
        if leftChr < rightChr:
            starts.append(columnNames[i+1])
            ends.append(columnNames[i])
    starts.insert(0,columnNames[0])
    ends.append(columnNames[-1])
    return(starts,ends)

#also finds chromosome cutoffs in datafiles where column names are slightly different
def getChromosomeCutoffsRedo(df):
    starts = list()
    ends = list()
    columnNames = df.columns
    a = [x.split(".") for x in columnNames]
    
    for i in range(len(a)-1):
        rightChr = int(a[i+1][0])
        leftChr = int(a[i][0])
        if leftChr < rightChr:
            starts.append(columnNames[i+1])
            ends.append(columnNames[i])
    starts.insert(0,columnNames[0])
    ends.append(columnNames[-1])
    return(starts,ends)

#imports data files
def dataframeListCreator(samplesList):
    methylationFiles = [x+"_cpgs_w1_s1_Meth.csv" for x in samplesList]
    CNVFiles = [x+"_cpgs_w1_s1_CNV.csv" for x in samplesList]

    #import the methylation files
    methDfList = list()
    
    for i in methylationFiles:
        curDf = pd.read_csv(i,index_col=0)
        methDfList.append(curDf)

    CNVDfList = list()
    #import the CNV files
    for i in CNVFiles:
        curDf = pd.read_csv(i,index_col=0)
        CNVDfList.append(curDf) 
        
    return(methDfList,CNVDfList)

#converts lists of dataframes to a stacked dataframe
def dataframeGenerator(inputList):
    df = pd.concat(inputList)
    return(df)

#chromosomes that are not amplified or deleted are labelled neutral
def grouper(x,chrList):
    if x not in chrList:
        x = 'Neutral'
    return(x)

#visualizes methylation data grouped by chromosome
def visualizerGroupedChr(samplesList,chrList,order,color,ylimits,title,output):
    methDfList,CNVDfList = dataframeListCreator(samplesList)
    methDf = dataframeGenerator(methDfList)
    
    starts,ends = getChromosomeCutoffs(methDfList[0])
    
    colNames = list()
    for j in range(len(starts)):
        colName = j+1
        methDf[colName] = methDf.loc[:,starts[j]:ends[j]].mean(axis=1)
        colNames.append(colName)
    methDfMean = methDf[colNames]
    methDfMean = methDfMean.stack()
    df = pd.DataFrame(methDfMean)
    df2 =df.reset_index()
    df2['group'] = df2['level_1'].apply(lambda x: grouper(x,chrList))

    colors = dict(zip(order,color))

    fig = plt.figure(figsize=[10,10])

    plt.suptitle(title)
    ax = fig.add_subplot(1,1,1)
    ax.hlines(0,-1,len(order),color='0.7',ls='--')
    ax = sns.boxplot(x = 'group', y=0, data = df2,palette=colors,order=order)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('% DNA methylation \n (Tumor - Normal)')
    ax.set_ylim(ylimits)
#    ax.set_yticks([-20,-10,0,10,20])
    
    plt.savefig(output)
    plt.show()

#runs Mann Whitney U test between chromosomes
def mannWhitneyUImplementer(samplesList,chrList):
    methDfList,CNVDfList = dataframeListCreator(samplesList)
    methDf = dataframeGenerator(methDfList)
    
    starts,ends = getChromosomeCutoffs(methDfList[0])
    
    colNames = list()
    for j in range(len(starts)):
        colName = j+1
        methDf[colName] = methDf.loc[:,starts[j]:ends[j]].mean(axis=1)
        colNames.append(colName)
    methDfMean = methDf[colNames]
    methDfMean = methDfMean.stack()
    df = pd.DataFrame(methDfMean)
    df2 =df.reset_index()
    df2['group'] = df2['level_1'].apply(lambda x: grouper(x,chrList))
    chrList.append('Neutral')   
    combos = list(combinations(chrList,2))
    outputDict = dict()
    for i in combos:
        stat,p = mannwhitneyu(df2[0][df2['group'] == i[0]],df2[0][df2['group'] == i[1]],alternative = 'two-sided')
        outputDict[i] = p
    df3 = pd.DataFrame(outputDict.items(),columns=['test','p-value'])
    return(df3)  

#returns sample size for Mann Whiteny U test
def mannWhitneyUImplementerCounts(samplesList,chrList):
    methDfList,CNVDfList = dataframeListCreator(samplesList)
    methDf = dataframeGenerator(methDfList)
    
    starts,ends = getChromosomeCutoffs(methDfList[0])
    
    colNames = list()
    for j in range(len(starts)):
        colName = j+1
        methDf[colName] = methDf.loc[:,starts[j]:ends[j]].mean(axis=1)
        colNames.append(colName)
    methDfMean = methDf[colNames]
    methDfMean = methDfMean.stack()
    df = pd.DataFrame(methDfMean)
    df2 =df.reset_index()
    df2['group'] = df2['level_1'].apply(lambda x: grouper(x,chrList))
    chrList.append('Neutral')   
    combos = list(combinations(chrList,2))
    outputDict = dict()
    for i in combos:
        x = df2[0][df2['group'] == i[0]]
        y = df2[0][df2['group'] == i[1]]
        outputDict[i] = [x,y]
        stat,p = mannwhitneyu(df2[0][df2['group'] == i[0]],df2[0][df2['group'] == i[1]],alternative = 'two-sided')
        outputDict[i] = p
    df3 = pd.DataFrame(outputDict.items(),columns=['test','counts'])
    return(df3)            

#GBM
os.chdir('GBM_w1_s1_final')
samplesGBMNew = ["MGH105A","MGH105B","MGH105C","MGH105D","MGH121", "MGH124","MGH129","MGH115","MGH122","MGH211"]

for sample in samplesGBM:
    visualizerGroupedChr([sample],[7,10],[7,'Neutral',10],['r','0.7', 'g'],[-0.25,0.25],'Average methylation by chromosome \n Among intergenic regions \n %s' %sample,"Methylation_By_Chromosome_%s_7_10.pdf" %sample)

visualizerGroupedChr(['MGH121'],[7,10],[7,'Neutral',10],['r','0.7', 'g'],[-0.25,0.25],'Average methylation by chromosome \n MGH121',"Methylation_By_Chromosome_MGH121_7_10.pdf")

dfList = []
for i in samplesGBMNew:
    print(i)
    df = mannWhitneyUImplementer([i],[7,10])
    df = df.set_index('test')
    dfList.append(df)
df = pd.concat(dfList,axis=1)
df.columns = samplesGBM
df.to_csv('Methylation_By_Chromosome_GBM.csv')

df = mannWhitneyUImplementer(['MGH121'],[7,10])
df = df.set_index('test')
df.to_csv('Methylation_By_Chromosome_MGH121_7_10.csv')

dfList = []
for i in samplesGBMNew:
    print(i)    
    df = mannWhitneyUImplementerCounts([i],[7,10])
    df = df.set_index('test')
    dfList.append(df)
df = pd.concat(dfList,axis=1)
df.columns = samplesGBM
df.to_csv('Methylation_By_Chromosome_GBM_counts.csv')

#IDH
os.chdir('IDH_w1_s1_final')
samplesIDHNoDD = ['IDH_MGH142_P1', 'IDH_MGH142_P2', 'MGH45', 'MGH64', 'MGH107', 'MGH135']
visualizerGroupedChr(['MGH107'],[7,8],[7,8,'Neutral'],['r','r','0.7'],[-0.25,0.25],'Average methylation by chromosome \n MGH107',"Methylation_By_Chromosome_MGH107_7_8.pdf")
visualizerGroupedChr(['MGH135'],[1,19],['Neutral',1,19],['0.7','g','g'],[-0.25,0.25],'Average methylation by chromosome \n MGH135',"Methylation_By_Chromosome_MGH135_1_19.pdf")
visualizerGroupedChr(['IDH_MGH142_P1'],[1,19],['Neutral',1,19],['0.7','g','g'],[-0.25,0.25],'Average methylation by chromosome \n MGH142_P1',"Methylation_By_Chromosome_MGH142_P1_1_19.pdf")
visualizerGroupedChr(['IDH_MGH142_P2'],[1,19],['Neutral',1,19],['0.7','g','g'],[-0.25,0.25],'Average methylation by chromosome \n MGH142_P2',"Methylation_By_Chromosome_MGH142_P2_1_19.pdf")

visualizerGroupedChr(['IDH_MGH142'],[1,19],['Neutral',1,19],['0.7','g','g'],[-0.25,0.25],'Average methylation by chromosome \n MGH142',"Methylation_By_Chromosome_MGH142_1_19.pdf")


df = mannWhitneyUImplementer(['MGH107'],[7,8])
df = df.set_index('test')
df.to_csv('Methylation_By_Chromosome_MGH107_7_8.csv')

df = mannWhitneyUImplementer(['MGH135'],[1,19])
df = df.set_index('test')
df.to_csv('Methylation_By_Chromosome_MGH135_1_19.csv')

df = mannWhitneyUImplementer(['IDH_MGH142_P1'],[1,19])
df = df.set_index('test')
df.to_csv('Methylation_By_Chromosome_IDH_MGH142_P1_1_19.csv')

df = mannWhitneyUImplementer(['IDH_MGH142_P2'],[1,19])
df = df.set_index('test')
df.to_csv('Methylation_By_Chromosome_IDH_MGH142_P2_1_19.csv')

df = mannWhitneyUImplementer(['IDH_MGH142'],[1,19])
df = df.set_index('test')
df.to_csv('Methylation_By_Chromosome_IDH_MGH142_1_19.csv')

