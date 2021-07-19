import pandas as pd
import os

#transform the input files to match the file formats given by
#https://github.com/broadinstitute/inferCNV/wiki/File-Definitions

samplesGBM = ["105A","105B","105C","105D",124,129,122,211,115,121]

#Raw Counts Matrix for Genes x Cells
#convert the raw scRNA counts files into .csv files
normalFile = '/GBM_normal.min2000.rsem.counts.noZeros.noMT.noRP.txt'
normalDf = pd.read_csv(normalFile,sep='\t',index_col=0)
normalCells = normalDf.columns.values

os.chdir('/Gene_counts')

for i in samplesGBM:
    if i == "105C":
        file = 'MGH105C_bis.rsem.counts.txt'
    else:
        file = 'MGH%s.rsem.counts.txt' %i
    #import the file and get the names of the tumor cells for the sample annotation file
    df = pd.read_csv(file,sep='\t',index_col=0)
    tumorCells = df.columns.values
    
    #merge on the normal data
    #if a cell is in the normal datafile, keep the data from the normal cell
    df = df.join(normalDf,how='inner',lsuffix='_left',rsuffix='_right')
    for j in df.columns.values:
        if j.split('_')[-1] == 'left':
            df = df.drop(j,axis=1)
        if j.split('_')[-1] == 'right':
            newName = j[:-6]
            df.rename(columns={j:newName}, inplace=True)
    
    df.to_csv('/inferCNV_MGH'+str(i)+'/MGH'+str(i)+'.rsem.counts_inferCNV.txt',sep='\t')
    
    #Sample annotation file
    cellNames = list(tumorCells)
    cellNames.extend(list(normalCells))
    
    #because cells can be in the tumor and normal data files, we need to 
    #deduplicate
    cellNames = list(set(cellNames))
    #if a cell is in the normal data file, it counts as a normal (even if it
    #is also in the tumor data file)
    cellTypes = ['normal' if x in normalCells else 'tumor' for x in cellNames]
    annotationDf = pd.DataFrame([cellNames,cellTypes])
    annotationDf = annotationDf.T
    annotationDf.to_csv('/inferCNV_MGH'+str(i)+'/MGH'+str(i)+'_annotation_file_inferCNV.txt',header=False,index=False,sep='\t')
    
    #gene ordering file
    genePosFile = '/MGH'+str(i)+'_gene_pos.txt'
    genePosDf = pd.read_csv(genePosFile,sep='\t',index_col=0)
        
    #delete chromosome 0
    genePosDf = genePosDf[genePosDf['chromosome_name'] > 0]
    genePosDf = genePosDf.drop('ensembl_gene_id',axis=1)
    genePosDf.to_csv('/inferCNV_MGH'+str(i)+'/MGH'+str(i)+'_gene_pos_inferCNV.txt',header=False,index=False,sep='\t')

