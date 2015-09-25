import os, sys

file_parts = {1:1986, 2:2206, 3:1864, 4:1973, 5:1707, 6:1622, 7:1473, 8:1359, 9:1013, 10:1141, 11:1186, 
             12:1230, 13:998, 14:823, 15:674, 16:657, 17:658, 18:702, 19:519, 20:481, 21:335, 22:273,
             'X':983, 'Y':24, 'M':1}


with open('/gpfs/group/torkamani/bhuvan/wellderly/coverage/all_samples.txt','r') as inp:
    lines = inp.readlines()

Samples = []
for line in lines:
    Samples.append(line.replace('\n',''))
    
datad = '/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageFiles'
efile = open('missing_coverages.txt','w',1)
for sample in Samples:
    subdir = os.path.join(datad,sample)
    print sample
    for chrm in range(1,23)+['X','Y','M']:
        print chrm
        sys.stdout.flush()
        numfiles = file_parts[chrm]
        counter  = 0
        
        while counter < numfiles:

            file_to_check = os.path.join(subdir,sample+'_'+str(chrm)+'_part_'+str(counter)+'.pickle')

            if not os.access(file_to_check,os.F_OK):
                efile.write(sample+'\t'+str(chrm)+'\t'+str(counter)+'\n')
                break
            counter += 1


