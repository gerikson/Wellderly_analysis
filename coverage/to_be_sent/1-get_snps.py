import sys,os
import cPickle
import gzip

	
def parse_line(line):
	
	line = line.split('\t')
	pos  = int(line[1])
	dlen = len(line[3])
	return pos, dlen

chrm = sys.argv[1]  ## Just the number or letter
ROOT = sys.argv[2]  ## Folder with the VCF files
OUTD = sys.argv[3]  ## Output Folder
print chrm
	
vcfFile = 'wellderly_inova.chr'+str(chrm)+'.vcf.gz'  #VCF File Name
vcfPath = os.path.join(ROOT,vcfFile)		     #VCF File Path

vcfHand = gzip.open(vcfPath)
snpLine = vcfHand.readline()

while not snpLine.startswith('#CHROM'):
	snpLine = vcfHand.readline()
		
snpLine = vcfHand.readline()


totalSNPs = 1000000
count     = 0
Filepart  = 0

while snpLine:	

	if count == 0:
		chrmDict = {}
		outFile  = 'snps_chrm_'+str(chrm)+'_part_'+str(Filepart)+'.pickle'   #Output File Name
		outPath  = os.path.join(OUTD,outFile)                 	        #Output File Path

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
		cPickle.dump(chrmDict,out)

	tln = 0
	for i in chrmDict:
		tln += len(chrmDict[i])

	Filepart += 1
	count = 0	



