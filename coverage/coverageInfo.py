import os, sys
import cPickle
import multiprocessing as mp
import pandas as pd
import functools
from pandas import DataFrame
from calculateCoverage import Coverage
import time

def mp_coverage_process(sampleName,sampleMap,snpsDict,sortedSnps):
	
		print sampleName
		sys.stdout.flush()	

		covrFile = sampleMap[sampleName]
		cinst    = Coverage(snpsDict,sortedSnps,covrFile)

		return sampleName, cinst.main()

class CoverageMain(object):

	def __init__(self,chrm,Filepart,datadir,outdir,snpsFilePath,SamplesAllFilePath,SamplesMedianFilePath):
		
		## snpsfile - This is the file for each chromosome with all the snps in it. It is a dictionary - {6:[1,3],10:[1],23:[1,5],....} (Pos:[lengths of variants at that pos])
		## 	      Produced by get_snps.py
		## SamplesFilePath - this file contains a list of all the DID sample Ids for samples from the white population, one on each line 

		self.chrm = 'chr'+chrm
		self.part = Filepart    ## Needed to generate a unique output file name
		self.ddir = datadir     ## Folder where all the pickled coverage files are
		self.odir = outdir	## Folder to write all the output files

		with open(snpsFilePath,'rb') as inp1:
			self.snps = cPickle.load(inp1)
		
		self.snpsPos = sorted(list(self.snps))

		with open(SamplesAllFilePath) as inp2:
			lines = inp2.readlines()
		samplesA  = []
		for line in lines:
			samplesA.append(line.replace('\n',''))
	
		self.samA = samplesA
		self.maps = self.get_files_for_chrm() 
		self.samA = self.filter_samples(samA)  ## Updating samples by removing the samples that are missing the coverages for the chromosome in question

		with open(SamplesMedianFilePath) as inp3:
			lines = inp3.readlines()
		samplesM  = []
		for line in lines:
			samplesA.append(line.replace('\n',''))
			
		self.samM = samples
		self.samM = self.filter_samples(samM)   ## Updating samples by removing the samples that are missing the coverages for the chromosome in question

	
	def main(self):
		
		sampleMap   = self.maps
		snpsDict    = self.snps
		MainDF      = DataFrame()
		snpsList    = self.snpsPos 
		samplesAll  = self.samA
		samplesUsed = self.samM

		## code for multiprocessing
		
		pool  = mp.Pool(processes=4)
		for res in pool.imap_unordered(functools.partial(mp_coverage_process,sampleMap=sampleMap,snpsDict=snpsDict,sortedSnps=snpsList),samplesAll,chunksize=50):
			sampleName, (coverages, snpIndices) = res
			DF_temp = DataFrame({sampleName:coverages},index = snpIndices)
			MainDF  = pd.concat([MainDF,DF_temp],axis = 1)
			
		"""
		for sampleName in samplesUsed:
			
			print sampleName
			sys.stdout.flush()	
			covrFile = sampleMap[sampleName]
			coverages, snpIndices = self.coverage_process(covrFile,snpsDict)
			DF_temp = DataFrame({sampleName:coverages},index = snpIndices)
			MainDF  = pd.concat([MainDF,DF_temp],axis = 1)
		
		"""		
		
		subDF   = MainDF.filter(items=samplesUsed)
		medians = subDF.median(axis=1)
		medFile = os.path.join(self.odir,'median_coverages_for_chrm_'+self.chrm+'_part'+self.part+'.csv')
		medians.to_csv(medFile)

		outf = os.path.join(self.odir,'coverages_for_chrm_'+self.chrm+'_part'+self.part+'.csv')
		MainDF.to_csv(outf)

		return 'Done'
	
	def get_files_for_chrm(self):
		
		chrm    = self.chrm
		samples = self.samA
		datadir = self.ddir

		mappings = {}
		for sample in samples:
			
			fname = sample+'_'+chrm+'.pickle'
			fpath = os.path.join(datadir,sample,fname)
			mappings[sample] = fpath
		
		return mappings

	def filter_samples(self,samplesAll):
		
		sampleMap   = self.maps

		samplesUsed = []
		for sampleName in samplesAll:
			covrFile = sampleMap[sampleName]
			if os.access(covrFile,os.F_OK):
				samplesUsed.append(sampleName)
			else:
				print 'Missing Sample = '+sampleName
				sys.stdout.flush()

		return samplesUsed		
		
	def coverage_process(self,coverageFile,snpsDict):
	
		cinst    = Coverage(snpsDict,coverageFile)
		return cinst.main()
		

if __name__ == '__main__':

	chrm  = sys.argv[1]
	fpart = sys.argv[2]
	ddir  = sys.argv[3]
	snpsP = sys.argv[4]
	
	#chrm  = '22'
	#fpart = '0'
	#ddir  = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageFiles'
	#snpsP = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/SnpFiles/snps_chrm_22_part_0.pickle'

	outd   = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo'
	samAll = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/'
	white  = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/White_0.95.txt'	
	
	t1 = time.time()
	cMain = CoverageMain(chrm,fpart,ddir,outd,snpsP,white)
	cMain.main()
	t2 = time.time()

	print t2-t1
	


