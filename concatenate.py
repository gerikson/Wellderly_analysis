import os, sys

file_parts = {'1':20, '2':23, '3':19, '4':20, '5':18, '6':17, '7':15, '8':14, '9':11, '10':12, '11':12,
             '12':13, '13':10, '14':9, '15':7, '16':7, '17':7, '18':8, '19':6, '20':5, '21':4, '22':3,
             'X':10, 'M':1}

dataDir = '/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/Medians'
chrm = sys.argv[1]

part = 0
outfile = '/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/MediansCompiled/medians_chrm_'+str(chrm)+'.csv'
while part <  file_parts[chrm]:
   infile = '/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/Medians/median_coverages_for_chrm_'+str(chrm)+'_part'+str(part)+'.csv'
   os.system('cat '+infile+' >> '+outfile)
   part += 1