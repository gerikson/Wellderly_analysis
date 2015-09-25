import os

for chrm in  range(1,23)+['X','Y','M']:

	with open('get_snps.job','w') as job:
		job.write('cd /gpfs/group/torkamani/bhuvan/wellderly/coverage\n')
		job.write('python get_snps.py '+str(chrm)+'\n')
	
	os.system('qsub get_snps.job')	


