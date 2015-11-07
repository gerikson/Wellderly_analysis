"""
Remove 24mer and ALL clustered variants

"""
import os, sys, gzip, datetime

class Overlap(object):
 
        def __init__(self,arr1Start,arr2Start,arr1End,arr2End):
 
        ## Two sets of start and end sites. Each start and end define the
        ## various islands
 
                self.a1s = arr1Start
                self.a2s = arr2Start
                self.a1e = arr1End
                self.a2e = arr2End


        def main(self):
 
                a1s = self.a1s
                a2s = self.a2s
                a1e = self.a1e
                a2e = self.a2e
 
                overlap_array = self.overlapCounts(a1s,a2s,a1e,a2e)
                return overlap_array
 

       	
       	#Overlaping one begin and end with and array of overlaps (a2s and a2e)
        def overlapCounts(self,a1s,a2s,a1e,a2e):
       
 
                #AllOverlaps = []
                #i = 0
                #while i < len(a1s):
                       
                totalOverlap = 0
                #beg1 = a1s[i]
                #end1 = a1e[i]
                beg1 = int(a1s)
                end1 = int(a1e)


                indx = self.insert_pos(a2s,a2s,beg1)
                if indx == "ERROR":
                    return False

                try:
                        beg2 = a2s[indx]
                        end2 = a2e[indx]
                except IndexError:
                       
                        indx = indx-1
                        beg2 = a2s[indx]
                        end2 = a2e[indx]
               
                overlap = self.overlapNumber(beg1,end1,beg2,end2)
                totalOverlap += overlap

                flagF, indxF = True, indx
                flagB, indxB = True, indx

                while flagF:

                        overlap = False
                        indxF += 1
                        try:
                                beg2 = a2s[indxF]
                                end2 = a2e[indxF]
                        except IndexError:
                                flagF = False
                                break

                        overlap = self.overlapNumber(beg1,end1,beg2,end2)
                        if overlap:
                                return True
                                #totalOverlap += overlap
                        else:
                                flagF = False

                while flagB:

                        overlap = False
                        indxB -= 1
                       
                        if indxB < 0:
                                break
                        try:
                                beg2 = a2s[indxB]
                                end2 = a2e[indxB]
                        except IndexError:
                                flagB = False
                                break
                       
                        overlap = self.overlapNumber(beg1,end1,beg2,end2)
                        if overlap:
                                return True
                                #totalOverlap += overlap
                        else:
                                flagB = False
                       
                       
 
                        #i += 1
 
                #print (AllOverlaps)
                return False   
 
        def insert_pos(self,arrI,arr,value):
 
        ## arrI - to get the final index
        ## arr  - the array to work on, same as arrI
       
            if len(arr) == 1:
            	if value > arr[0]:
            		return arrI.index(arr[0])+1
            	else:
            		return arrI.index(arr[0])
            else:
            	arr_mid = len(arr)/2
            	arrLow  = arr[:arr_mid]
            	arrHigh = arr[arr_mid:]
                #print arrLow
            	#try:
            	if value>arrLow[-1]:
            		out = self.insert_pos(arrI,arrHigh,value)
            	else:
            		out = self.insert_pos(arrI,arrLow,value)

                return out


               
        def overlapPercent(self,s1,e1,s2,e2):
 
        ## The function which calculates the overlap between two isalnds
        ##
        ## seq1  s1---------------------e1
        ## seq2             s2---------------------e2
        ##
        ## seq1  s1----------------------------e1
        ## aeq2            s2--------------e2
        ##
        ## seq1           s1--------------------e1
        ## aeq2    s2----------------e2
        ##
        ## seq1            s1-----------e1
        ## aeq2    s2-----------------------------e2
        ##
                
				if (s1<=s2) and (e1<e2) and (e1>=s2):
					return True
				elif (s1<=s2) and (e1>=e2):
					return True
				elif (s1>s2) and (s1<=e2) and (e1>=e2):
					return True
				elif (s1>s2) and (e1<e2):
					return True
				else:
					return False

 				'''
                seqL1 = float(e1-s1+1)
                seqL2 = float(e2-s2+1)
 
                if (s1<=s2) and (e1<e2) and (e1>=s2):
                        overlap = float(e1-s2+1)
                        ovrPerc = min((overlap/seqL1),(overlap/seqL2))*100.0
                elif (s1<=s2) and (e1>=e2):
                        overlap = seqL2
                        ovrPerc = (overlap/seqL1)*100.0
                elif (s1>s2) and (s1<=e2) and (e1>=e2):
                        overlap = float(e2-s1+1)
                        ovrPerc = min((overlap/seqL1),(overlap/seqL2))*100.0           
                elif (s1>s2) and (e1<e2):
                        overlap = seqL1
                        ovrPerc = (overlap/seqL2)*100.0
                else:
                        ovrPerc = 0.0
 
                return ovrPerc 
                '''
 



        def overlapNumber(self,s1,e1,s2,e2):
 
			seqL1 = float(e1-s1+1)
			seqL2 = float(e2-s2+1)
				

			if (s1<=s2) and (e1<e2) and (e1>=s2):
			        overlap = float(e1-s2+1)
			elif (s1<=s2) and (e1>=e2):
			        overlap = seqL2
			elif (s1>s2) and (s1<=e2) and (e1>=e2):
			        overlap = float(e2-s1+1)
			elif (s1>s2) and (e1<e2):
			        overlap = seqL1
			else:
			        overlap = 0.0
			return overlap   
 


def main(chrom):

    ch = "1"

    kmer_file="/gpfs/home/gerikson/wellderly/filter_data/24kmer/byChrom/unique_24mer" + str(chrom) + ".txt"
    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/vcf_rareVariant."+chrom+".vcf.gz"
    
    output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/24mer_ALLcluster_filtres/vcf_rareVariant."+chrom+".vcf.gz"
    clustered_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_ALL_clustered/ALL_clustered."+chrom+".txt.gz"
    counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/24mer_ALLcluster_filtres/counter.txt"

    repeatMask_begin = []
    repeatMask_end = []


    print 'Prepping Kmers...'
    count_rep = 0

    print 'Prepping Kmers...'
    with open(kmer_file) as f:
        for line in f:
            items = line.split('\t')
            repeatMask_begin.append(int(items[1]))
            repeatMask_end.append(int(items[2]))
    
    
    print "Prepping clustered file"
    clustered_dict = {}
    with open(clustered_file) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                clustered_dict[items[1]] = items[3]



    print 'Calculating...'

    counter = 0
    good_lines = 0


    repeatMask_counter = 0
    block = ""

    f = gzip.open(input_filename)
    o = gzip.open(output_filename, 'w')
    #filter_block = ""

    for line in f:


        if line[:1] == "#":
            print "header"
            o.write(line)
            continue

        counter += 1

        if counter%100 == 0:
        	o.write(block)
        	#filt.write(filter_block)
        	block = ""
        	#filter_block = ""
        	o.flush()

        if counter%10000 == 0:
        	print datetime.datetime.now().time()
        	print "total lines " + str(counter)
        	print "Good lines " + str(good_lines)
        	sys.stdout.flush()


        tp_line = line.strip().split("\t")

        begin = int(tp_line[1])
        end = int(tp_line[1]) + 1

        #Verify if this is one of the clustered variants
        try:
            var = clustered_dict[str(begin)]
            #If this has a different start
            if var != tp_line[3]:
                #Check k_mers
                repOverlap = Overlap(begin, repeatMask_begin, end, repeatMask_end)
                overlapRep = repOverlap.main()
                if overlapRep:
                    block = block + line
                    good_lines += 1
                    continue
        except:
            #Check k_mers
            repOverlap = Overlap(begin, repeatMask_begin, end, repeatMask_end)
            overlapRep = repOverlap.main()
            if overlapRep:
                block = block + line
                good_lines += 1
                continue
        
    		

    #filt.write(filter_block)
    o.write(block)

    #print "Total g counter " + str(repeatMask_counter)

    print "Total lines " + str(counter)
    print "Total good lines " + str(good_lines)
    c = open(counter_file, "a")
    c.write(chrom+"\t"+str(counter)+"\t"+str(good_lines)+"\n")
    c.close()
    f.close()
    o.close()
    #filt.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)