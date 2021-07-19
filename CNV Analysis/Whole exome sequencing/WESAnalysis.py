import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

#importt the exome files
os.chdir('/Users/lkluegel/Documents/GBM/exome')
exomeSample = ['105A','105B','105C','122','115','124V2']

sample = exomeSample[0]
file = 'MGH%sExome.csv' %sample
df = pd.read_csv(file,index_col=0)

#only keep relevant columns
df = df[['Chromosome','Start.bp','End.bp','total_copy_ratio']]

'''in the datafiles, each chromosome starts at 0bp. for plotting purposes,
the previous chromosomes end #bp is added to the next to get a nice axis
the next two functions do this'''
def endLocRebaser(row):
    curChr = row['Chromosome']
    if curChr > 1:
        prevChrs = [x for x in range(1,int(curChr))]
        newBase = 0
        for i in prevChrs:
            newBase = newBase + df[df['Chromosome'] == i].iloc[-1,2]
        curValue = row['End.bp']
        output = curValue + newBase
    else:
        output = row['End.bp']
    return(output)
    
def startLocRebaser(row):
    curChr = row['Chromosome']
    if curChr > 1:
        prevChr = curChr-1
        newBase = df[df['Chromosome'] == prevChr].iloc[-1,2]
        curValue = row['Start.bp']
        output = curValue + newBase
    else:
        output = row['Start.bp']
    return(output)
    



df['End.bp'] = df.apply(lambda x: endLocRebaser(x),axis=1)
df['Start.bp'] = df.apply(lambda x: startLocRebaser(x),axis=1)

#plot the data
os.chdir('/Users/lkluegel/Documents/GBM/exome')
for sample in exomeSample:
    file = 'MGH%sExome.csv' %sample
    df = pd.read_csv(file,index_col=0)
    df = df[['Chromosome','Start.bp','End.bp','corrected_total_cn']]
    df['End.bp'] = df.apply(lambda x: endLocRebaser(x),axis=1)
    df['Start.bp'] = df.apply(lambda x: startLocRebaser(x),axis=1)
    df['corrected_total_cn'] = df['corrected_total_cn'].apply(lambda x: x if x>0.2 else 0.2)
    fig, ax = plt.subplots()
    for i in range(df.shape[0]):
        row = df.iloc[i,:]
        curValue = np.log2(row['corrected_total_cn']/2)
        if curValue < -0.3:
            c = 'b'
        elif curValue > 0.3:
            c = 'r'
        else:
            c = 'k'
        ax.hlines(y=np.log2(row['corrected_total_cn']/2),xmin=row['Start.bp'],xmax=row['End.bp'],color ='k')
        if row['End.bp'] - row['Start.bp'] <4000000:
            ax.plot(row['Start.bp'],np.log2(row['corrected_total_cn']/2),marker='.',color='k',markersize=1)
    xticks = list(df.groupby('Chromosome').min()['Start.bp'])
    xticks.append(max(df['End.bp']))
    ax.vlines(x=xticks,ymin=[min(np.log2(df['corrected_total_cn']/2)) for x in range(23)],ymax=[max(np.log2(df['corrected_total_cn']/2)) for x in range(23)],linewidth=0.1)
    xlabels = [x if x%2 != 0 else "" for x in range(1,23)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('log2(copy number)')
    plt.title('CNV based on WES \n MGH%s' %sample)
    plt.savefig('/Users/lkluegel/Documents/GBM/exome/WES_MGH%s.pdf' %sample)
    plt.show()
