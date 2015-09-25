import bz2
import sys, os, cPickle

class EofF(object):
    
    def __init__(self):
        pass
    def __eq__(self, other):
        return False
    def __ne__(self, other):
        return True
    def __lt__(self, other):
        return False
    def __gt__(self, other):
        return True
    def __le__(self, other):
        return False
    def __ge__(self,other):
        return False


class CoverageParser(object):

	def __init__(self,coverageFilePath):

		self.covr = bz2.BZ2File(coverageFilePath,'rb')
		self.pInd, self.cInd = self.get_coverage_file_header()

	def __iter__(self):
		
		return self
	
	def next(self):

		line = self.get_coverage_file_line()
		if line == (None,None):
			raise StopIteration
		else:
			return line
	
	def get_coverage_file_header(self):

		fhand = self.covr
		line  = fhand.readline()
		
		while not line.startswith('>offset'):
			line = fhand.readline()
		
		line = line.split('\t')
		cInd = line.index('uniqueSequenceCoverage')
		pInd = line.index('>offset')

		return pInd, cInd   ## Indices for the position and coverage in coverage file
	
	def get_coverage_file_line(self):

		fhand = self.covr
		line  = fhand.readline()
		
		if line == '':
			return EofF,EofF()
		else:
			line = line.split('\t')
			posn = int(line[self.pInd])+1
			covr = int(line[self.cInd])
			return posn, covr
	
	def get_coverage_file_block(self,blockSize):
		
		c = []
		p = []
		i = 0
		while i< blockSize:
			posn, covr = self.get_coverage_file_line()
			c.append(covr)
			p.append(posn)
			i += 1
		return p,c


class Coverage(object):

	def __init__(self,snpsDict,sortedSnps,coverageFilePath):

		self.snps  = snpsDict			## Dictionary with all the snps and lengths
		self.slist = sortedSnps			## Sorted list of all snp positions
		scratchDr  = os.environ['PBSTMPDIR']

		CFP, CFN   = os.path.split(coverageFilePath)
		self.NCFP  = os.path.join(scratchDr,CFN)              #New coverage file path
		os.system('cp '+coverageFilePath+' '+self.NCFP)
		
		with open(self.NCFP,'rb') as inp:
			self.covr = cPickle.load(inp)

	def main(self):
		
		coverages  = []
		snpIndices = []

		snpDict = self.snps
		snpList = self.slist
			
		for snpPos in snpList:
			for snpLen in snpDict[snpPos]:
				
				snpInd = str(snpPos)+':'+str(snpLen)
				snpCov = self.get_snp_coverage(snpPos,snpLen)
				
				coverages.append(snpCov)
				snpIndices.append(snpInd)	
		
		os.remove(self.NCFP)
		return coverages, snpIndices

	def get_snp_coverage(self,pos,length):

		covD  = self.covr
		bases = []
		i = 0
		cover = 0
		while i<length:
			try:
				cover += covD[pos+i]
			except IndexError:
				print 'Coverage Missing at Position ' + str(pos)
			i += 1
			
		aCovr = float(cover)/length
		return aCovr


if __name__ == '__main__':

    t = EofF()
    print 100000000==t
     
    





				
				


			











			
			
		


	
