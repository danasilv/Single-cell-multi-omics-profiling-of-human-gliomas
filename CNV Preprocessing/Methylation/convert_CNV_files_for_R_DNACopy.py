import pandas as pd
import os

'''this script gets the CNV data into the right shape for DNACopy'''

os.chdir('./Full_genome')

#import the single cell datafiles of choice
samplesGBM = ["105A","105B","105C","105D",124,129,122,211,115,121]
dfList = list()
for i in samplesGBM:
    file = 'GBM_MGH%s_cpgs_w10_s1_CNV.csv' %i
    df = pd.read_csv(file,index_col=0)
    df['origin'] = ['MGH'+str(i) for x in df.index]
    dfList.append(df)
    
dfNew = pd.concat(dfList)

#average by sample 
dfAvg = dfNew.groupby('origin').mean()

#get the dataframe into the shape required
dfAvg = dfAvg.T
cols = [x for x in dfAvg.columns.values]

dfAvg.reset_index(level=0, inplace=True)

#DNACopy needs a location column, which here is just the index
cols = ['location'] + cols
dfAvg.columns = cols

new = dfAvg.location.str.split(pat='.',expand=True)
dfAvg['Chromosome'] = new[0]
dfAvg['Position'] = new[1]
dfAvg = dfAvg.drop(columns = 'location')

cols = ['Chromosome','Position']+['MGH'+str(x) for x in samplesGBM]
dfAvg = dfAvg[cols]
dfAvg.columns = cols

#here, the data is subset for chromosome 7, this application dependent
dfAvgChr7 = dfAvg[dfAvg['Chromosome'] == '7']
dfAvgChr7 = dfAvgChr7.reset_index()

dfAvgChr7['Position'] = [x+1 for x in dfAvgChr7.index]

dfAvg = dfAvg[['Chromosome','Position']+['MGH'+str(x) for x in samplesGBM]]

dfAvg.to_csv('Chr7_CNV_for_R_DNAcopy_w10_s1.csv',index=False)

