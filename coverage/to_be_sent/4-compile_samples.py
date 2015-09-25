import pandas as pd
from pandas import DataFrame 
import os, sys, cPickle

chrm  = sys.argv[1]  
part  = sys.argv[2]  ## This is the partition number of the variants
datad = sys.argv[3]  ## Directory where all the files containing the coverages are
sampA = sys.argv[4]  ## A file containing the sample IDs that you want to compile together (One ID on each line)
sampM = sys.argv[5]  ## A file containing a subset of the samples for which you want to calculate the median coverage per variant (we are using the whites)
outd  = sys.argv[6]  ## Directory to store the results 
             ## This directory should have 3 sub folders - All, Whites, Medians

with open(sampA,'r') as inp:
    lines = inp.readlines()
samplesAll = []
for line in lines:
    samplesAll.append(line.replace('\n',''))

with open(sampM,'r') as inp:
    lines = inp.readlines()
samplesMed = []
for line in lines:
    samplesMed.append(line.replace('\n',''))

MainDF = DataFrame()
for sample in samplesAll:

    subDirectory = os.path.join(datad,sample)
    covrFileName = sample+'_'+chrm+'_part_'+str(part)+'.pickle'
    covrFilePath = os.path.join(subDirectory,covrFileName)

    try:    
        with open(covrFilePath,'rb') as inp:
            data = cPickle.load(inp)
    except IOError,e:
        print e

    coverages  = data['coverages']
    varIndices = data['indices']

    DF_temp = DataFrame({sample:coverages},index = varIndices)
    MainDF  = pd.concat([MainDF,DF_temp],axis = 1)


subDF   = MainDF.filter(items=samplesMed)
medians = subDF.median(axis=1)
medFile = os.path.join(outd,'Medians','median_coverages_for_chrm_'+chrm+'_part'+part+'.csv')
medians.to_csv(medFile)

outf1 = os.path.join(outd,'Whites','coverages_for_chrm_'+chrm+'_part'+part+'_whites_only.csv')
subDF.to_csv(outf1)
outf2 = os.path.join(outd,'All','coverages_for_chrm_'+chrm+'_part'+part+'_all_samples.csv')
MainDF.to_csv(outf2)





    



