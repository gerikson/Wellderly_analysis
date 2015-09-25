import sys,os
import cPickle
import gzip

ROOT = '/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes'  ## Folder with the VCF files
OUTD = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/SnpFiles'				       ## Output Folder	

def parse_line(line):
	
	line = line.split('\t')
	pos  = int(line[1])
	dlen = len(line[3])
	return pos, dlen

#for chrm in range(1,23)+['X','Y','M']:
#for chrm in ['M']:
chrm = sys.argv[1]
print chrm
	
vcfFile = 'wellderly_inova.chr'+str(chrm)+'.vcf.gz'  #VCF File Name
vcfPath = os.path.join(ROOT,vcfFile)		     #VCF File Path

vcfHand = gzip.open(vcfPath)
snpLine = vcfHand.readline()

while not snpLine.startswith('#CHROM'):
	snpLine = vcfHand.readline()
		
snpLine = vcfHand.readline()


totalSNPs = 10000
count     = 0
Filepart  = 0

posFile  = 'snps_chrm_pos_'+str(chrm)+'.txt'	  #Output file to save positions of variants
posFileP = os.path.join(OUTD,posFile)
posFileH = open(posFileP,'wb')

while snpLine:	

	if count == 0:
		chrmDict = {}
		outFile  = 'snps_chrm_'+str(chrm)+'_part_'+str(Filepart)+'.pickle' 	  #Output File Name
		outPath  = os.path.join(OUTD,outFile)                 	       	   	  #Output File Path

	while (count<totalSNPs) and (snpLine):
		
		pos,dlen = parse_line(snpLine)	
		if pos in chrmDict:
			if dlen in chrmDict[pos]:
				pass
			else:
				chrmDict[pos].append(dlen)
				chrmDict[pos] = sorted(chrmDict[pos])
		else:
			chrmDict[pos] = [dlen]
		snpLine = vcfHand.readline()
		count += 1

	with open(outPath,'w') as out:
		cPickle.dump(chrmDict,out,2)

	beg = str(min(chrmDict))               ## Position of the first variant
	end = max(chrmDict)	  	       ## Position of the last variatnt
	end = str(end+(max(chrmDict[end])-1))  ## End position of the largest variant at the last start position 

	posFileH.write(('\t').join([str(Filepart),beg,end])+'\n')

	Filepart += 1
	count = 0	

posFileH.close()

