import os
import pandas as pd

'''path to CBS.py file'''
os.chdir(pathToCBS.py)
import CBS

'''file to datafiles'''
os.chdir(pathToDataFiles)

file1 = 'EGFR_MGH105A_CNV_Chr7.csv'
file2 = 'EGFR_MGH105B_CNV_Chr7.csv'
file3 = 'EGFR_MGH105C_CNV_Chr7.csv'
file4 = 'EGFR_MGH105D_CNV_Chr7.csv'
file5 = 'EGFR_MGH115_CNV_Chr7.csv'
file6 = 'EGFR_MGH121_CNV_Chr7.csv'
file7 = 'EGFR_MGH122_CNV_Chr7.csv'
file8 = 'EGFR_MGH124_CNV_Chr7.csv'
file9 = 'EGFR_MGH129_CNV_Chr7.csv'
file10 = 'EGFR_MGH211_CNV_Chr7.csv'


df1 = pd.read_csv(file1, index_col = 0)
df2 = pd.read_csv(file2, index_col = 0)
df3 = pd.read_csv(file3, index_col = 0)
df4 = pd.read_csv(file4, index_col = 0)
df5 = pd.read_csv(file5, index_col = 0)
df6 = pd.read_csv(file6, index_col = 0)
df7 = pd.read_csv(file7, index_col = 0)
df8 = pd.read_csv(file8, index_col = 0)
df9 = pd.read_csv(file9, index_col = 0)
df10 = pd.read_csv(file10, index_col = 0)


frames = [df1,df2,df3,df4,df5,df6,df7,df8,df8,df9,df10]

df = pd.concat(frames)

samples = ['GBM_MGH105A','GBM_MGH105B','GBM_MGH105C','GBM_MGH105D','GBM_MGH115','GBM_MGH121','GBM_MGH122','GBM_MGH124','GBM_MGH129','GBM_MGH211']

df['new_col'] = df.index
df['origin'] = df['new_col'].str[4:11].map(lambda x: x.rstrip('_'))
del df['new_col']

output1 = CBS.outer(df,0.999,10000,True,5272355760701257514)

CBS.CBSVisualizerGrouped(df,output1,'CNVs based on circular binary segmentation \n Chromosome 7','output.pdf','log')

