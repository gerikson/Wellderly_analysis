"""
Extract header from the annotation file

"""
import os, sys, gzip, datetime

def main(input_filename):

    counter =0

    line_buffer = ""
    header = ""
    #head = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/HEADER.wellderly_inova_annotations_part3")
    head = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/HEADER.batch_WHspocvQya_smallvars_vypdIYlUQa_FINAL")
    for line in head:
        counter += 1
        if counter%100 == 0:

            outfile.write(line_buffer)
            line_buffer = ""
            outfile.flush()


        if counter%100000 == 0:
            print str(counter)


        templine=line.strip().split("\t")
        print templine[1:7]
        
        for index, i in enumerate(templine):
            if i == 'Gene':
                print "Gene " + str(index)
            if i == 'Coding_Impact':
                print "Index coding impact " + str(index)
            if i == '1000genomes_EUR':
                print "Index 1000genomes_EUR " + str(index)
        

  
    
if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)