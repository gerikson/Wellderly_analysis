"""
Split annotation file by chromosome

"""
import os, sys, gzip, datetime

def main(input_filename):

    inp = gzip.open(input_filename)
    outfilename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/byChrom/"
    
    chrom="chr1"
    name=outfilename+chrom+".txt.gz"
    outfile=gzip.open(name,"a")
    counter =0

    line_buffer = ""
    header = ""
    #head = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/HEADER.wellderly_inova_annotations_part3")
    
    for line in inp:
    #for line in head:
        counter += 1
        if counter%100 == 0:

            outfile.write(line_buffer)
            line_buffer = ""
            outfile.flush()


        if counter%100000 == 0:
            print str(counter)
            print datetime.datetime.now().time()
            sys.stdout.flush()


        templine=line.strip().split("\t")
        '''
        print templine[1:7]
        
        for index, i in enumerate(templine):
            if i == 'Coding_Impact':
                print "Index coding impact " + str(index)
            if i == '1000genomes_EUR':
                print "Index 1000genomes_EUR " + str(index)
            if i == 'Splice_Site_Pred':
                print "Splice_Site_Pred " +str(index)
            if i == 'Location':
                print "Location " +str(index)
        '''
        if counter != 1:

            try:
                if chrom == templine[0]:
                    #Extracting only the first 6 columns plus Gene, Coding Impact and 1000Genomes_EUR
                    #line_buffer = line_buffer + "\t".join(templine[1:7])+ "\t"+templine[8]+"\t"+templine[12]+"\t"+templine[37] +"\n"
                    #For the 3rd file:
                    line_buffer = line_buffer + "\t".join(templine[:6])+ "\t"+templine[20]+"\t"+templine[24]+"\t"+templine[49] +"\t"+templine[22]+"\t"+templine[73]+"\n"
                else:

                    outfile.write(line_buffer)
                    line_buffer = ""
                    outfile.flush()
                    outfile.close()
                    #print "New Chrom "+ templine[0]

                    chrom=templine[0]
                    name=outfilename+chrom+".txt.gz"
                    outfile=gzip.open(name,"a")
                    #line_buffer = "\t".join(templine[1:7])+"\t"+templine[8]+"\t"+templine[12]+"\t"+templine[37] +"\n"     
                    #For the 3rd file:
                    line_buffer = line_buffer + "\t".join(templine[:6])+ "\t"+templine[20]+"\t"+templine[24]+"\t"+templine[49] +"\t"+templine[22]+"\t"+templine[73]+"\n"   
            except:
                print line
        
    outfile.write(line_buffer)
    line_buffer = ""
    outfile.flush()
    outfile.close()
    inp.close()
  
    
if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)