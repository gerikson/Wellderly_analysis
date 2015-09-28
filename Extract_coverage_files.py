import os
import cPickle

Root = '/gpfs/group/stsi/data/projects/wellderly'

output_wg="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/total_cov_wg_and_exome.txt"
output_wg="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/whiteOnly_wg_and_exome.txt"
whitesfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/White_0.95.txt"

outp = open(output_wg, 'w')

whites_id = {}
for i in wh:
	whites_id[i] = "Y"


DirecDict = {}
for batchNum in ['01','02','03','04','05','06','07','08','09','10','11']:
	batchName   = 'batch'+batchNum
	batchFolder = os.path.join(Root,batchName)
	
	print batchName

	for dirs, subdirs, files in os.walk(batchFolder):
		#print files
		for f in files:
			if 'summary' in f:
				#print f
				t_dir = dirs.split("/")
				#print t_dir[7] + "\t" + f
				final_line = t_dir[7] + "\t" + f
				opfile = open(f)

				for line in opfile:
					if 'Genome coverage Gross mapping yield (Gb)' in line:
						final_line = final_line + "\t" + line.strip()
					if 'Exome coverage  Exome fraction where weightSumSequenceCoverage >= 10x' in line:
						final_line = final_line + "\t" + line.strip()

				try:
					#
					whites_id[t_dir[7]]

				#print dirs

		'''
		for subdir in subdirs:
			if 'ASM' == subdir:
				

				
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
		'''

'''
with open('coverageFilePaths.pickle','wb') as outp:
	cPickle.dump(DirecDict,outp)
'''