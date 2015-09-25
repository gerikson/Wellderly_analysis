import os, sys

file_parts = {'1':1986, '2':2206, '3':1864, '4':1973, '5':1707, '6':1622, '7':1473, '8':1359, '9':1013, '10':1141, '11':1186, 
             '12':1230, '13':998, '14':823, '15':674, '16':657, '17':658, '18':702, '19':519, '20':481, '21':335, '22':273,
             'X':983, 'Y':24, 'M':1}

dataDir = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo/Medians'
chrm = sys.argv[1]

part = 0
outfile = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo/MediansCompiled/medians_chrm_'+str(chrm)+'.csv'
while part <  file_parts[chrm]:
   infile = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo/Medians/median_coverages_for_chrm_'+str(chrm)+'_part'+str(part)+'.csv'
   os.system('cat '+infile+' >> '+outfile)
   part += 1
