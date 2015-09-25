import os, sys
import cPickle
from calculateCoverage import CoverageParser

fpath   = sys.argv[1]
sample  = sys.argv[2]
chrm    = sys.argv[3]  ## eg 1, 2 ,3, X, M
outd    = sys.argv[4]
datadir = sys.argv[5]

subd = os.path.join(outd,sample)
os.system('mkdir '+subd)
print sample, chrm

class FileNotFound(Exception):
    pass

def get_snpsDict(chrm,filepart):

    snpsFileName = 'snps_chrm_'+chrm+'_part_'+str(part)+'.pickle'
    snpsFilePath = os.path.join(datadir,snpsFileName)
    if os.access(snpsFilePath,os.F_OK):
        print chrm,filepart
        with open(snpsFilePath,'rb') as inp:
            varDict = cPickle.load(inp)         ## This dictionary is one of the file parts of the main VCF file containing 10,000 variants
        varList = sorted(list(varDict))             ## This is the sorted position of the chromosome where the Variants are
        return varDict,varList
    else:
        raise FileNotFound('End of fileparts reached')

def get_covr(pos,length,posL,covL):

    ## This function takes in a position, the length of the variant at that position,
    ## a list of positions and a list of coverages at those positions and returns the avg coverage.
    ## A position that doesn't have an entry in the data file is considered to have 0 coverage
    
    i = 0
    tCov = 0
    tLen = 0
    while i < length:
        try:
            tCov = tCov+covL[posL.index(pos+i)]
            tLen += 1
            i += 1
        except ValueError:
            i += 1
    if tLen != 0: cov = tCov/float(tLen)
    else: cov = 0
    return cov

def cov_for_missing_pos(pos,posL,covL):
    
    # This function calculates the coverage for a variant for which the 1st bp is missing
    i = 1
    tCov = 0
    tLen = 0
    while i<35:
        try:
            tCov = tCov+covL[posL.index(pos+i)]
            tLen += 1
            try:
                tCov = tCov+covL[posL.index(pos-i)]
                tLen += 1
                i += 1

            except ValueError:
                i+= 1
        except ValueError:
            try:
                tCov = tCov+covL[posL.index(pos-i)]
                tLen += 1
                i += 1
            except ValueError:
                i+= 1
    
    if tLen != 0: cov = tCov/float(tLen)
    else: cov = 0

    return cov


## Starting with File part 0
part = 0
varDict,varList = get_snpsDict(chrm,part)

d = {'indices':[],'coverages':[]}
bufferSize = 200
bufferPos  = [None]*bufferSize
bufferCov  = [None]*bufferSize

#from pandas import DataFrame as DF
cpIn = CoverageParser(fpath)

while True:
    pos,covr = cpIn.get_coverage_file_line()
    bufferPos.pop(0)
    bufferPos.append(pos)
    bufferCov.pop(0)
    bufferCov.append(covr)
    
    
    if varList == []:
        
        ## We have reached the end of this file part
        ## Writing out the local coverage dictionary and updating it
        
        outf = os.path.join(subd,sample+'_'+chrm+'_part_'+str(part)+'.pickle')
        with open(outf,'wb') as outp:
            cPickle.dump(d,outp,2)
        d = {'indices':[],'coverages':[]}
        
        ## Updating Snps Dictionary
        part = part+1
        try:
            varDict,varList = get_snpsDict(chrm,part)  
        except FileNotFound,e:
            print e
            break
    
    if  varList[0] > pos:
        ## This means we can skip this position as there isn't any new variant here
        pass

    elif  varList[0] == pos:
        
        #print bufferPos
        
        varPos    = pos 
        variants  = varDict[varPos]
        blockSize = max(variants)-1
        posns, covrs = cpIn.get_coverage_file_block(blockSize)
            
        positions = [varPos] + posns    
        coverages = [covr] + covrs
        ##Saving Coverages in the dictionary
        for varLen in variants:
            varID    = str(varPos)+'_'+str(varLen)
            varCov   = get_covr(varPos,varLen,positions,coverages)
            d['indices'].append(varID)
            d['coverages'].append(varCov)
        
        ## Updating bufer
        bufferPos = bufferPos[blockSize:]
        bufferPos.extend(posns)
        bufferCov = bufferCov[blockSize:]
        bufferCov.extend(covrs)
        
        ## Removing the topmost Variant Position 
        varList.pop(0)
            
    else: # varList[0] < pos:
        
        ## This means that we have gone further and we'll have to check in the buffer
        varPos    = varList[0]
        variants  = varDict[varPos]
        
        try:
            bufferInd = bufferPos.index(varPos)
            len_remaining_entries = (bufferSize-bufferInd)
            
            if max(variants) <= len_remaining_entries:
                blockSize = 0
            else:
                blockSize = max(variants)-len_remaining_entries
    
            posns, covrs = cpIn.get_coverage_file_block(blockSize)      
            coverages = bufferCov[bufferInd:] + covrs
            positions = bufferPos[bufferInd:] + posns

            ##Saving Coverages in the dictionary
            for varLen in variants:
                varID    = str(varPos)+'_'+str(varLen)
                varCov   = get_covr(varPos,varLen,positions,coverages)  
                d['indices'].append(varID)
                d['coverages'].append(varCov)
            
            ## Updating bufer
            bufferPos = bufferPos[blockSize:]
            bufferPos.extend(posns)
            bufferCov = bufferCov[blockSize:]
            bufferCov.extend(covrs)
            
            ## Removing the topmost Variant Position 
            varList.pop(0)

        except ValueError:
            # This is in case the coverage file doesn't have coerage for this base
            # We take the average of a read length (70 bp) around the position
            
            rightEnd  = varPos+35
            blockSize = rightEnd - bufferPos[-1]
            if blockSize>0:
                
                ## Updating bufer if right side 35 bp are not present in the buffer already
                posns, covrs = cpIn.get_coverage_file_block(blockSize)
                bufferPos = bufferPos[blockSize:]
                bufferPos.extend(posns)
                bufferCov = bufferCov[blockSize:]
                bufferCov.extend(covrs)
            else:
                pass
            
            varCov = cov_for_missing_pos(varPos,bufferPos,bufferCov)                
            ##Saving Coverages in the dictionary
            for varLen in variants:
                varID    = str(varPos)+'_'+str(varLen)
                varCov   = get_covr(varPos,varLen,bufferPos,bufferCov)  
                d['indices'].append(varID)
                d['coverages'].append(varCov)
    
            ## Removing the topmost Variant Position 
            varList.pop(0)

    

