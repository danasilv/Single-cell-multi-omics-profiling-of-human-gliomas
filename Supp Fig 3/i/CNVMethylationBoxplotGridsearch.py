import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.stats import mannwhitneyu

#threshold CNV status
def CNVStatus(x,threshold):
    if x > threshold:
        status = 1
    elif x < -threshold:
        status = -1
    else:
        status = 0
    return(status)

#create the downsampled files
os.chdir('/Users/lkluegel/Documents/GBM/CNV_data/GBM_w1_s1_final')
samplesGBM = ["MGH105A","MGH105B","MGH105C","MGH105D",'MGH121_P1','MGH121_P2', 'MGH121_P3', 'MGH121_P4', 'MGH124','MGH129','MGH115','MGH122','MGH211']
windows = [1,3,5,10]
thresholdList = [0.3,0.5,0.7,0.9]

#runs the methylation and CNV downsampling scripts
def dataframeDownsampledCreator(samplesList, windows):
    for sample in samplesList:
        for window in windows:
            print(sample)
            print(window)
            os.system('python MethylationDownsampled.py GBM_w1_s1_final %s_cpgs_w1_s1_methCpGs.txt %s_cpgs_w1_s1.txt %s_cpgs_w1_s1_Meth_window_%s.csv GBM %s' %(sample, sample, sample, window, window))
            os.system('python CNVDownsampled.py GBM_w1_s1_final  %s_cpgs_w1_s1.txt %s_cpgs_w1_s1_CNV_window_%s.csv GBM %s' %(sample, sample, window, window))
    
dataframeDownsampledCreator(samplesGBM, windows)

#loads the methylation and CNV files and gets them in the right shape for downstream analysis
def dfLoader(samplesList, window):
    methylationFiles = [x+"_cpgs_w1_s1_Meth_window_%s.csv" %window for x in samplesList]
    CNVFiles = [x+"_cpgs_w1_s1_CNV_window_%s.csv" %window for x in samplesList]
    
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

#stacks the columns of interest to get the data in the right shape for boxplots
def dataframeTransformerNew(methDfList,CNVDfList,threshold):
    Y = pd.concat(methDfList)
    
    Y = Y.unstack()

    X = pd.concat(CNVDfList)
    
    X = X.unstack()

    df = X.to_frame().merge(Y.to_frame(),left_index=True,right_index=True)
    df.columns = ['CNV','Methylation']

    #remove methylation outliers
    YQ1 = df['Methylation'].quantile(0.25)
    YQ3 = df['Methylation'].quantile(0.75)
    IQR = YQ3-YQ1
    
    df = df[df['Methylation']>YQ1-1.5*IQR]
    df = df[df['Methylation']<YQ3+1.5*IQR]
    #limit CNV
    df = df[df['CNV']<=6]

    df['CNVStatus'] = df['CNV'].apply(lambda x:CNVStatus(x,threshold))
    return(df)

#draws the boxplot grids
def visualizerNew(samplesList,windowList,thresholdList,title,output):
    numWindows = len(windowList)
    numThresholds = len(thresholdList)
    fig = plt.figure(figsize=[10,10])

    plt.suptitle(title)
    ax = fig.add_subplot(1,1,1)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    my_pal = {1: "r", 0: "0.7", -1:"g"}    
    num=1
    for i in range(numWindows):
        methDfList,CNVDfList = dfLoader(samplesList, windowList[i])
        for j in range(numThresholds):
            ax = fig.add_subplot(numWindows,numThresholds,num)
            df = dataframeTransformerNew(methDfList,CNVDfList,thresholdList[j])
            ax = sns.boxplot(data = df, x = 'CNVStatus', y = 'Methylation', palette = my_pal, order=[1,0,-1], showfliers=False)
            ax.set_xticklabels(labels=['Ampl.','Neutral','Del.'])
            ax.set_ylim([-0.65,0.65])
            ax.set_title('win = %s, thresh = %s' %(windowList[i],thresholdList[j]))
            ax.set_xlabel("")
            ax.set_ylabel("")

            num += 1
    fig.tight_layout(rect=[0, 0.03, 1, 0.9])

    plt.savefig(output)
    plt.show()

#run the visualizer
visualizerNew(samplesGBM,windows,thresholdList,"Methylation by copy number aberration \n GBM samples","Methylation_by_CNV_gridsearch_GBM.pdf")

#get the Mann-Whitney test scores
tests = ['Amp. vs Neu.','Amp vs Del','Del vs Neu']
columns = pd.MultiIndex.from_product([thresholdList, tests], names=['Threshold', 'Test'])
outputDf = pd.DataFrame(index = windows, columns = columns)

for i in range(len(windows)):
    methDfList,CNVDfList = dfLoader(samplesGBM, windows[i])
    for j in range(len(thresholdList)):
        df = dataframeTransformerNew(methDfList,CNVDfList,thresholdList[j])
        x = df['Methylation'][df['CNVStatus']==1]
        y = df['Methylation'][df['CNVStatus']==0]
        z = df['Methylation'][df['CNVStatus']==-1]
        score1 = mannwhitneyu(x,y)
        score2 = mannwhitneyu(x,z)
        score3 = mannwhitneyu(y,z)
        outputDf.loc[windows[i],(thresholdList[j],tests[0])] = score1[1]
        outputDf.loc[windows[i],(thresholdList[j],tests[1])] = score2[1]
        outputDf.loc[windows[i],(thresholdList[j],tests[2])] = score3[1]

outputDf.to_csv('Methylation_by_CNV_gridsearch_GBM_MannTest_scores.csv')

#Do the same things for IDH
os.chdir('IDH_w1_s1_final')
samplesIDHNoDD = ['IDH_MGH142_P1', 'IDH_MGH142_P2', 'MGH45', 'MGH64', 'MGH107', 'MGH135']

def dataframeDownsampledCreator(samplesList, windows):
    for sample in samplesList:
        for window in windows:
            print(sample)
            print(window)
            os.system('python MethylationDownsampled.py IDH_w1_s1_final %s_cpgs_w1_s1_methCpGs.txt %s_cpgs_w1_s1.txt %s_cpgs_w1_s1_Meth_window_%s.csv IDH %s' %(sample, sample, sample, window, window))
            os.system('python CNVDownsampled.py IDH_w1_s1_final  %s_cpgs_w1_s1.txt %s_cpgs_w1_s1_CNV_window_%s.csv IDH %s' %(sample, sample, window, window))
    
dataframeDownsampledCreator(samplesIDHNoDD, windows)

visualizerNew(samplesIDHNoDD,windows,thresholdList,"Methylation by copy number aberration \n IDH samples without double digest","Methylation_by_CNV_gridsearch_IDH.pdf")

tests = ['Amp. vs Neu.','Amp vs Del','Del vs Neu']
columns = pd.MultiIndex.from_product([thresholdList, tests], names=['Threshold', 'Test'])
outputDf = pd.DataFrame(index = windows, columns = columns)

for i in range(len(windows)):
    methDfList,CNVDfList = dfLoader(samplesIDHNoDD, windows[i])
    for j in range(len(thresholdList)):
        df = dataframeTransformerNew(methDfList,CNVDfList,thresholdList[j])
        x = df['Methylation'][df['CNVStatus']==1]
        y = df['Methylation'][df['CNVStatus']==0]
        z = df['Methylation'][df['CNVStatus']==-1]
        score1 = mannwhitneyu(x,y)
        score2 = mannwhitneyu(x,z)
        score3 = mannwhitneyu(y,z)
        outputDf.loc[windows[i],(thresholdList[j],tests[0])] = score1[1]
        outputDf.loc[windows[i],(thresholdList[j],tests[1])] = score2[1]
        outputDf.loc[windows[i],(thresholdList[j],tests[2])] = score3[1]

outputDf.to_csv('Methylation_by_CNV_gridsearch_IDH_no_double_digest_MannTest_scores.csv')
