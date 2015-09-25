import os

file_parts = {1:1986, 2:2206, 3:1864, 4:1973, 5:1707, 6:1622, 7:1473, 8:1359, 9:1013, 10:1141, 11:1186, 
             12:1230, 13:998, 14:823, 15:674, 16:657, 17:658, 18:702, 19:519, 20:481, 21:335, 22:273,
             'X':983, 'Y':24, 'M':1}

for chrm in [10]: #file_parts:
    with open('concat.job','w') as o:
        o.write('cd /gpfs/group/torkamani/bhuvan/wellderly/coverage\n')
        o.write('python concatenate.py '+str(chrm)+'\n')
    os.system('qsub concat.job')
