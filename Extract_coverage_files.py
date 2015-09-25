import os
import cPickle

Root = '/gpfs/group/stsi/data/projects/wellderly'

DirecDict = {}
for batchNum in ['01','02','03','04','05','06','07','08','09','10','11']:
	batchName   = 'batch'+batchNum
	batchFolder = os.path.join(Root,batchName)
	
	print batchName

	for dirs, subdirs, files in os.walk(batchFolder):
		for subdir in subdirs:
			if 'REF' in subdir:
				
				coverageDir = os.path.join(dirs,'REF')
				temp = coverageDir.split(os.sep)

				didID = temp[8].split('-')[0]
				asmID = temp[9]
					
				print didID	
				
				DirecDict[didID] = {}
				for chrm in range(1,23)+['X','Y','M']:

					coverageFile = 'coverageRefScore-chr'+str(chrm)+'-'+asmID+'.tsv.bz2'
					coveragePath = os.path.join(coverageDir,coverageFile)

					DirecDict[didID]['chr'+str(chrm)] = coveragePath



with open('coverageFilePaths.pickle','wb') as outp:
	cPickle.dump(DirecDict,outp)
