import os, cPickle

with open('/gpfs/group/torkamani/bhuvan/wellderly/coverage/coverageFilePaths.pickle','rb') as inp:
    FileMaps = cPickle.load(inp)
#with open('/gpfs/group/torkamani/bhuvan/wellderly/coverage/all_samples.txt','r') as inp:
#    lines = inp.readlines()

## Running for only missed samples

with open('missing_coverages.txt','r') as inp:
    lines = inp.readlines()

MissedSamples = {}
for line in lines:
    line = line.split('\t')
    samp = line[0]
    chrm = line[1]
    
    if samp in MissedSamples:
        MissedSamples[samp].append(chrm)
    else:
        MissedSamples[samp] = [chrm]  

Samples = []
for line in lines:
    Samples.append(line.replace('\n',''))
    
outd = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageFiles'
efile = open('missing_chrms_B.txt','w')

#for sample in ['GS000026838']:
for sample in MissedSamples:
    
    print sample        
    #for chrm in [21,14]:
    for chrm in MissedSamples[sample]:
        print chrm

        chrm  = str(chrm)
        fpath = FileMaps[sample]['chr'+chrm]
        datad = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/SnpFiles'

        if os.access(fpath,os.F_OK):

            with open('dictionaryB.job','w') as outp:
                
                outp.write('#!/bin/csh\n')
                outp.write('module load python\n')
                outp.write('cd /gpfs/group/torkamani/bhuvan/wellderly/coverage\n')
                outp.write('python make_coverage_dicts.py '+fpath+' '+sample+' '+str(chrm)+' '+outd+' '+datad+'\n')

            os.system('qsub dictionaryB.job')
        else:
            efile.write(sample+'\t'+chrm+'\t'+fpath+'\n')
        
efile.close()


