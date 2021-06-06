import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

samplesGBM = ["105A","105B","105C","105D",124,129,122,211,115,121]

'''get the gene expression output generated by inferCNV
these files come from inferCNVBySample.R'''
EGFRDfList = list()
for i in samplesGBM:
    os.chdir('/Users/lkluegel/Documents/GBM/scRNA/inferCNV_MGH'+str(i))
    
    #import the non-smoothed gene expression file
    file = 'preliminary_modified_expression.csv'
    df = pd.read_csv(file,
                     index_col=0)

    #import the annotation file
    annotationFile = 'MGH'+str(i)+'_annotation_file_inferCNV.txt'
    annotationDf = pd.read_csv(annotationFile,
                               sep='\t',
                               index_col=0,
                               header=None)
    annotationDf.columns=['cell_type']
    
    #join it to the CNV file and remove the normal cells
    df = df.T
    df = df.join(annotationDf)
    df = df[df['cell_type'] != 'normal']
    
    #now drop the new column and transpose back
    df = df.drop('cell_type',axis=1)
    df = df.T
    
    #import the gene position file
    genePosFile = 'MGH'+str(i)+'_gene_pos_inferCNV.txt'
    genePosDf = pd.read_csv(genePosFile,
                            sep='\t',
                            index_col=0,
                            header=None)
    genePosDf.columns = ['chromosome','start','end']
    
    #join the data frames
    df = df.join(genePosDf,how='inner')
    
    #only keep data of interest (on chr7) and find the EGFR locus
    chr7 = df[df['chromosome'] == 7]
    EGFRStart = df.loc['EGFR','start']
    EGFREnd = df.loc['EGFR','end']

    #find the genes that are 25mb in either direction of the EGFR locus and 
    #remove the rest
    EGFR = chr7[chr7['start'] >= EGFRStart - 25000000]
    EGFR = EGFR[EGFR['end'] <= EGFREnd + 25000000]
    EGFR = EGFR.sort_values(by='start',axis=0)
    
    #drop the unnecessary columns
    EGFR = EGFR.drop(['chromosome','start','end'],axis=1)

    EGFRDfList.append(EGFR)

#take the first dataframe
EGFRDf = EGFRDfList[0]

'''iteratively join the other dataframes. This has to be joined not concatenated
because different dfs can contain different genes'''
for i in EGFRDfList[1:]:
    EGFRDf = EGFRDf.join(i,
                         how='inner')

EGFRDf = EGFRDf.T
EGFRDf['origin'] = [x.split('_')[0] for x in EGFRDf.index]
EGFRDf.to_csv('/Users/lkluegel/Documents/GBM/scRNA/Chr7_EGFR_25mb_inferCNV_All_GBM.csv')

#combine by origin to get the output data for DNAcopy in R
EGFRDfFlat = EGFRDf.groupby('origin').mean()
EGFRDfFlat = EGFRDfFlat.T
EGFRDfFlat['Chromosome'] = [7 for x in EGFRDfFlat.index]
EGFRDfFlat['Position'] = [x+1 for x in range(EGFRDfFlat.shape[0])]

#turn the EGFRDfFlat dataframe to log2 form for DNAcopy
for column in EGFRDfFlat.columns.values:
    if column[:3] == 'MGH':
        EGFRDfFlat[column] = EGFRDfFlat[column].apply(np.log2)
        
#this file went into DNAcopy
EGFRDfFlat.to_csv('/Users/lkluegel/Documents/GBM/scRNA/Chr7_EGFR_25mb_inferCNV_for_DNAcopy_All_GBM.csv')

'''open the DNAcopy output files with the breakpoints
these files come from CNV_CBS_EGFR_25mb.R'''
os.chdir('/Users/lkluegel/Documents/GBM/scRNA/DNAcopy_output')

#create the dictionary that will house the cutoff points
outputDict = dict()

#open each cutoff file
for i in samplesGBM:
    file = 'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH'+str(i)+'.csv'
    df = pd.read_csv(file,
                     index_col=0)
    
    #get the cutoffs
    cutoffs = df.loc[:,'loc.start']
    
    #shift them one to the left because R doesn't use 0-based indexing
    cutoffs = [x-1 for x in cutoffs]
    
    #add the last location
    cutoffs.append(max(df.loc[:,'loc.end'])-1)
    outputDict['MGH'+str(i)] = cutoffs
      
#now convert to the input data to log form and visualize
EGFRDfLog = EGFRDf.copy(deep=True)
for column in EGFRDfLog.columns.values[:-1]:
    EGFRDfLog[column] = EGFRDfLog[column].apply(np.log2)

'''used to visualize grouped or segmented data
copied the function here from CBS.py to change some of the parameters
it also returns a returnDict object that contains the 
CNV cutoff points'''
def CBSVisualizerGrouped2(dataframe,cutoffsDict,title,output):
    
    lastColumn = dataframe.columns.values[-1]
    subclonesByFreq = list(dataframe[lastColumn].value_counts().index)

    fig = plt.figure(figsize=[5,2*len(subclonesByFreq)+1])
    plt.suptitle(title,x=0.55)
    
    #the first set of axes is to provide the location for the main axis labels
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Chromosome location')
    ax.set_ylabel('Log2(Inferred CNV)')

    #hide the axes
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    returnDict = dict()
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
        axNew.vlines('EGFR',ymin=-3.5,ymax=6.5,color='k')
        axNew.plot(dfSummary['Mean CNV'],'-',color='0.4',linewidth=0.5)

        #only show the labels and ticks of the cutoff points to avoid clutter
        xtickLabels = ['']*cutoffs[-1]
        for i in cutoffs[1:-1]:
            dataLabel = dfSummary.index[i+1]
            xtickLabels[i] = dataLabel
        xtickLabels[:] = [x for x in xtickLabels if x != '' ]
        returnDict[j] = xtickLabels
        xticks = [x for x in xtickLabels if len(x) > 0]

        axNew.set_xticks(xticks)
        axNew.set_xticklabels([])
        #depending on the scale of the data, these cutoffs might need to be adjusted
        axNew.set_yticks([-1,0,1])
        axNew.set_ylim([-1.5,1.5])
        axNew.hlines(0,0,cutoffs[-1],'0.5',linewidth=0.1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.show()
    return(returnDict)
    
returnDict = CBSVisualizerGrouped2(EGFRDfLog,outputDict,'CNV based on scRNA using inferCNV','/Users/lkluegel/Documents/GBM/scRNA/Chr7_EGFR_25mb_scRNA_inferCNV_All_GBM_Log2_rescaled.pdf')
   