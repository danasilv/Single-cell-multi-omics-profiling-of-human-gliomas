'''This code is a simple implementation of the circular binary segmentation 
algorithm described by Olshen & Lucito (Circular binary segmentation for the 
analysis ofarray-based DNA copy number data)'''

import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
    
'''calculate the partial sums as described in the paper'''
def partialSums(dataframe,i,j,n):
    if i == 0:
        Si = dataframe.iloc[:,i].mean()
    else:
        Si = dataframe.iloc[:,:i+1].mean().sum()
    Sj = dataframe.iloc[:,:j+1].mean().sum()
    Sn = dataframe.mean().sum()
    return(Si,Sj,Sn)

'''calculate the Z values as described in the paper'''
def zGenerator(Si,Sj,Sn,i,j,n):
    Z1 = (1/(j-i)+1/(n-j+i))**(-1/2)
    Z2 = (Sj-Si)/(j-i)
    Z3 = (Sn-Sj+Si)/(n-j+i)
    Z = Z1*(Z2-Z3)
    return(Z)

'''adds new elements to the search queue'''
def searchQExtender(searchQueue,inputList,element):
    elementIndex = inputList.index(element)
    leftSearchRegion = [inputList[elementIndex-1],inputList[elementIndex]]
    rightSearchRegion = [inputList[elementIndex],inputList[elementIndex+1]]
    searchQueue.extend([leftSearchRegion,rightSearchRegion])
    return(searchQueue)

'''this implements a Monte Carlo simulation to find the cutoff point (ZCutoff) 
to find statistically significant differences between candidate ranges'''
def monteCarloZFinder(dataframe,sigLevel,numSimulations,seed):
    random.seed(seed)
    n = dataframe.shape[1]
    ZList = list()
    for k in range(numSimulations):
        i = random.randint(0,n-2)
        j = random.randint(i+1,n-1)
        Si, Sj, Sn = partialSums(dataframe,i,j,n)
        Z = zGenerator(Si,Sj,Sn,i,j,n)
        ZList.append(Z)
    ZCutoff = np.quantile(ZList,sigLevel)
    return(ZCutoff)

'''for a given region, determines if there is a significantly different region
Returns the most significantly different region  ''' 
def cutoffFinder(dataframe):
    n = dataframe.shape[1]
    maxZ = 0
    besti = 0
    bestj = 0
    for i in range(n-1):
        for j in range(i+1,n):
            Si, Sj, Sn = partialSums(dataframe,i,j,n)
            Z = zGenerator(Si,Sj,Sn,i,j,n)
            if abs(Z) > maxZ:
                maxZ = abs(Z)
                besti = i
                bestj = j    
    return(besti,bestj,maxZ)

'''takes a window of the original data from the search queue and runs 
cutoffFinder on that region. If there is a significantly different region
in that, it adds the boundaries of the regions to a list and generates new 
smaller windows that are added to the search queue'''
def CBS(dataframe,ZCutoff,searchQ,cutoffs):
    searchWindow = searchQ.pop()
    curDataframe = dataframe.iloc[:,searchWindow[0]:searchWindow[1]+1]
    besti, bestj, maxZ = cutoffFinder(curDataframe)
    #add recursion and add the known cutoffs to a list
    if maxZ > ZCutoff:
        bestiNew = 0
        bestjNew = 0
        if besti not in cutoffs:
            cutoffs.append(besti)
            bestiNew = 1
        if bestj not in cutoffs:
            cutoffs.append(bestj)
            bestjNew = 1
        if bestiNew + bestjNew > 0:
            cutoffs = sorted(cutoffs)
            if bestiNew == 1:
                searchQ = searchQExtender(searchQ,cutoffs,besti)
            if bestjNew == 1:
                searchQ = searchQExtender(searchQ,cutoffs,bestj)

            if bestiNew + bestjNew == 2:
                searchQ = set(tuple(x) for x in searchQ)
                searchQ = [list(x) for x in searchQ]
    return(searchQ,cutoffs)

'''runs all of the above functions while there are windows to be searched 
(i.e. there are values in the search queue) , return the list of cutoffs'''
def runCBS(dataframe,ZCutoff):
    n = dataframe.shape[1]
    cutoffs = [0,n]
    #instantiate the first search window
    searchQ = [[0,n]]
    while len(searchQ) > 0:
        searchQ,cutoffs = CBS(dataframe,ZCutoff,searchQ,cutoffs)
    
    if len(cutoffs) > 2:
        if cutoffs[1] == cutoffs[0] + 1:
            del cutoffs[1]
        if len(cutoffs) > 2:
            if cutoffs[-2] == cutoffs[-1] - 1:
                del cutoffs[-2]
    return(cutoffs)

'''wrapper function for runCBS to allow setting parameters for significance
levels, number of Monte Carlo Simulations, and random seed. The group parameter
allows data to be split out by subclone or sample. For that, the last column
has to contain sample/subclone name'''
def outer(dataframe,sigLevel,numSimulations,groups,seed):
    if groups == True:
        cutoffsDict = dict()
        lastColumn = dataframe.columns.values[-1]
        subclonesByFreq = list(dataframe[lastColumn].value_counts().index)
        for i in subclonesByFreq:
            curDf = dataframe[dataframe[lastColumn] == i].iloc[:,:-1]
            ZCutoff = monteCarloZFinder(curDf,sigLevel,numSimulations,seed)
            cutoffs = runCBS(curDf,ZCutoff)
            cutoffsDict[i] = cutoffs
        return(cutoffsDict)
    elif groups == False:
        ZCutoff = monteCarloZFinder(dataframe,sigLevel,numSimulations,seed)
        cutoffs = runCBS(dataframe,ZCutoff)
        return(cutoffs)
    

    
'''used to visualize grouped or segmented data from runCBS and outer functions'''
def CBSVisualizerGrouped(dataframe,cutoffsDict,title,output,method):
    
    lastColumn = dataframe.columns.values[-1]
    subclonesByFreq = list(dataframe[lastColumn].value_counts().index)

    fig = plt.figure(figsize=[5,2*len(subclonesByFreq)+1])
    plt.suptitle(title,x=0.55)
    
    #the first set of axes is to provide the location for the main axis labels
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Chromosome location')
    if method == 'log':
        ax.set_ylabel('Inferred CNV \n Log2(Tumor/Normal)')
    elif method == 'nolog':
        ax.set_ylabel('Inferred CNV \n Turmor/Normal')
    #hide the axes
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
    #make a subplot for each group
    for j in subclonesByFreq:
        subclone = subclonesByFreq.index(j)+1
        cutoffs = cutoffsDict[j]
        
        axNew = fig.add_subplot(len(subclonesByFreq),1,subclone)
        axNew.set_title(str(j),fontsize=12)
        curDf = dataframe[dataframe[lastColumn] == j].iloc[:,:-1]
        dfMeans = curDf.mean()
        region = list()
        for i in range(len(cutoffs)-1):
            localRegion = [i]*(cutoffs[i+1]-cutoffs[i])
            region.extend(localRegion)
        region.insert(0,0)
        dfSummary = pd.DataFrame(zip(dfMeans,region))
        dfSummary.index = dfMeans.index
        dfSummary.columns = ['Mean CNV','Region']
        dfSummary['Region mean'] = dfSummary['Mean CNV'].groupby(dfSummary['Region']).transform('mean')
        plt.plot()
        #plot the original individual data
        for i in range(curDf.shape[0]):
            plt.plot(curDf.iloc[i,:],'0.7',linewidth=0.1)
            
        #plot the average of by segment
        y = dfSummary['Mean CNV'].groupby(dfSummary['Region']).mean()
        xmin = list()
        xmax = list()
        for i in range(len(y)):
            xStart = cutoffs[i]+1
            xEnd = cutoffs[i+1]+1   
            xmin.append(xStart)
            xmax.append(xEnd)
        xmin[0] = 0
        xmax[-1] = cutoffs[-1]-1
        axNew.hlines(y,xmin,xmax,colors='r',zorder=10,linewidth=0.7)

        axNew.plot(dfSummary['Mean CNV'],'-',color='0.4',linewidth=0.5)
        
        #only show the labels and ticks of the cutoff points to avoid clutter
        xtickLabels = ['']*cutoffs[-1]
        for i in cutoffs[1:-1]:
            dataLabel = dfSummary.index[i+1]
            xtickLabels[i] = dataLabel
        xticks=[]
#        xticks = [x for x in xtickLabels if len(x) > 0]
#        xtickLabels = [x.split('.')[0] for x in xtickLabels if len(x)>0]

#        xtickLabels = [x.split('.')[1] for x in xtickLabels if len(x)>0]
        axNew.set_xticks(xticks)
#        axNew.set_xticklabels(xtickLabels)
        #depending on the scale of the data, these cutoffs might need to be adjusted
        if method == 'log':
            axNew.set_yticks([-3,0,3,6])
            axNew.set_ylim([-3.5,6.5])
#            axNew.hlines(0,0,cutoffs[-1],'0.5',linewidth=0.1)
        elif method == 'nolog':
             axNew.set_yticks([0,1,2])
             axNew.set_ylim([-0.5,2.5])   
             axNew.hlines(1,0,cutoffs[-1],'0.5',linewidth=0.1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.show()
 
