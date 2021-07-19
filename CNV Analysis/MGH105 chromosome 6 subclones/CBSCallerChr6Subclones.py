import pandas as pd
import os

os.chdir(pathToCBS.py)
import CBS
os.chdir(pathToDataFiles)

fileA = 'MGH105A_CNV_full_genome.csv'
fileB = 'MGH105B_CNV_full_genome.csv'
fileC = 'MGH105C_CNV_full_genome.csv'
fileD = 'MGH105D_CNV_full_genome.csv'

dfA = pd.read_csv(fileA, index_col = 0)
dfB = pd.read_csv(fileB, index_col = 0)
dfC = pd.read_csv(fileC, index_col = 0)
dfD = pd.read_csv(fileD, index_col = 0)

frames = [dfA,dfB,dfC,dfD]



df = pd.concat(frames)

subclonefile = 'GBM_CNV_MGH105_subclones_01212020.txt'
subcloneDf = pd.read_csv(subclonefile,sep='\t',index_col=0)

df = pd.merge(df,subcloneDf['subclone'],left_index = True, right_index = True)
subclones = set(df['subclone'])

chr6 = df[list(df.loc[:,'chr6.216':'chr6.249']) + ['subclone']]

output1 = CBS.outer(chr6,0.99,10000,True,5272355760701257514)

CBS.CBSVisualizerGrouped(chr6,output1,'CNVs based on circular binary segmentation \n By sub-clone: chromosome 6 MHG 105A/B/C/D','CNV_CBS_MGH105ABCD_Subclone.pdf','log')

