import os, sys
import cPickle
import multiprocessing as mp
import bz2
import pandas as pd
import functools
from pandas import DataFrame
from calculateCoverage import Coverage
import time

def coverage_process(sampleName,samples,snpsDict):
	
		print sampleName	
		#samples  = self.maps
		#snpsDict = self.snps
		covrFile = samples[sampleName]

		cinst = Coverage(snpsDict,covrFile)
		return sampleName, cinst.main()

class CoverageMain(object):

	def __init__(self,chrm,outdir,snpsFilePath,coverageDictFilePath):
		
		## snpsfile - This is the file for each chromosome with all the snps in it. It is a dictionary - {6:[1,3],10:[1],23:[1,5],....} (Pos:[lengths of variants at that pos])
		## coverageDictFilePath - This file contains a dictionary which maps every sample to the directory where the coverage is stored

		self.chrm = 'chr'+chrm
		self.odir = outdir
		self.maps = self.get_files_for_chrm(coverageDictFilePath) 	
		with open(snpsFilePath,'rb') as inp:
			self.snps = cPickle.load(inp)

	def main(self):
		
		samples  = self.maps
		snpsDict = self.snps
		MainDF   = DataFrame()

		errFile  = os.path.join(self.odir,'missin_coverage_file_'+self.chrm+'.txt')
		efile = open(errFile,'w')
		#temp = list(samples)
		
		samplesUsed = []
		for sampleName in samples:
			covrFile = samples[sampleName]
			if os.access(covrFile,os.F_OK):
				samplesUsed.append(sampleName)
			else:
				efile.write(sampleName+'\n')

		pool = mp.Pool()
		results = pool.map(functools.partial(coverage_process,samples=self.maps,snpsDict=self.snps),samplesUsed)
		
		"""
		for sampleName in samples:
			sys.stdout.flush()
			covrFile = samples[sampleName]
			
			if os.access(covrFile,os.F_OK):			
			
				coverages, snpIndices = self.coverage_process(snpsDict,covrFile)
				DF_temp = DataFrame({sampleName:coverages},index = snpIndices)
				MainDF  = pd.concat([MainDF,DF_temp],axis = 1)
			else:
				efile.write(sampleName+'\n')
		"""
		
		for res in results:

			sampleName, (coverages, snpIndices) = res
			DF_temp = DataFrame({sampleName:coverages},index = snpIndices)
			MainDF  = pd.concat([MainDF,DF_temp],axis = 1)

		efile.close()
		outf = os.path.join(self.odir,'coverages_for_chrm_'+chrm+'.csv')
		MainDF.to_csv(outf)

		return MainDF
				
		
	def get_files_for_chrm(self,cfile):
		
		chrm = self.chrm
		with open(cfile,'rb') as inp:
			data = cPickle.load(inp)
	
		mappings = {}
		for sample in data:
			mappings[sample] = data[sample][chrm]
	
		del data
		return mappings

		

if __name__ == '__main__':

	covrD = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/coverageFilePaths.pickle'
	chrm  = sys.argv[1]
	outd  = '/gpfs/group/torkamani/bhuvan/wellderly/coverage'
	snpsF = 'all_snps_chrm_'+chrm+'.pickle'
	snpsP = os.path.join(outd,snpsF)
		
	t1 = time.time()
	cMain = CoverageMain(chrm,outd,snpsP,covrD)
	cMain.main()
	t2 = time.time()

	print t2-t1
	


