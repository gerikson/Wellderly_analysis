"""
Count repeats_homopolymers etc by part

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
					try:
						if value>arrLow[-1]:
							out = self.insert_pos(arrI,arrHigh,value)
						else:
							out = self.insert_pos(arrI,arrLow,value)
					except:
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
 


def main(chrom, part):

	#Repeats
	#repeats_filename="/gpfs/home/gerikson/wellderly/filter_data/SimpleRepeats/byChrom/simpleRepeats."+chrom+".txt"
	#Micrositelites
	#repeats_filename="/gpfs/home/gerikson/wellderly/filter_data/Micrositelite/byChrom/micrositelite."+chrom+".txt"
	#Homopolymer
	#repeats_filename="/gpfs/home/gerikson/wellderly/filter_data/Homopolymer/byChrom/homopolymer."+chrom+".txt"
	#SegDup
	repeats_filename="/gpfs/home/gerikson/wellderly/filter_data/SegmentalDup/byChrom/SegmentalDup."+chrom+".txt"
	#RepeatMask
	#repeats_filename = "/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom/RepeatMasker."+chrom+".txt"

	repeats_begin = []
	repeats_end = []
	'''
	# Prep
	print 'Prepping repeats...'
	with open(repeats_filename) as f:
		for line in f:
			items = line.split('\t')
			repeats_begin.append(int(items[2]))
			repeats_end.append(int(items[3]))
	
	
	# Prep
	print 'Prepping microsite...'
	with open(repeats_filename) as f:
		for line in f:
			items = line.split('\t')
			repeats_begin.append(int(items[2]))
			repeats_end.append(int(items[3]))

	
	# Prep
	print 'Prepping homopolymer...'
	with open(repeats_filename) as f:
		for line in f:
			items = line.split('\t')
			repeats_begin.append(int(items[1]))
			repeats_end.append(int(items[2]))
	
	'''
	print 'Prepping SegDup...'
	with open(repeats_filename) as f:
		for line in f:
			items = line.split('\t')
			repeats_begin.append(int(items[2]))
			repeats_end.append(int(items[3]))
		
	'''
	print 'Prepping RepeatMask...'
	with open(repeats_filename) as f:
		for line in f:
			items = line.split('\t')
			repeats_begin.append(int(items[6]))
			repeats_end.append(int(items[7]))
	'''

	counter = 0
	counter_unique_no_pass = 0
	good_lines = 0
	
	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0


	AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/AF_parts/"+chrom+".part"+str(part)+".txt.gz")
	#filter_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/chr4_chr5_chr6_parts/repeat"+chrom+".part"+str(part)+".txt.gz", "w")
	#filter_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/chr4_chr5_chr6_parts/microsit"+chrom+".part"+str(part)+".txt.gz", "w")
	#filter_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/chr4_chr5_chr6_parts/homopolymer"+chrom+".part"+str(part)+".txt.gz", "w")
	filter_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/chr4_chr5_chr6_parts/SegDup"+chrom+".part"+str(part)+".txt.gz", "w")
	#filter_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/chr4_chr5_chr6_parts/RepeatMask"+chrom+".part"+str(part)+".txt.gz", "w")
	
	for line in AF_file:
		counter += 1
		tp_line = line.strip().split("\t")

		begin = int(tp_line[1])
		end = int(int(tp_line[1])) + int(len(tp_line[3]))
		
		#ADD repeats
		repeatOverlap = Overlap(begin, repeats_begin, end, repeats_end)
		overlapRep = repeatOverlap.main()
		if overlapRep:
			filter_file.write(line)
			counter_unique_no_pass += 1
			well_AF = float(tp_line[5])
			inova_AF = float(tp_line[6])
			#print str(well_AF)
			#print str(inova_AF)
			if well_AF > 0.0:
				if well_AF < 0.01:
				    well_001 += 1
				elif well_AF < 0.05:
				    well_005 += 1
				else:
				    well_005_plus += 1

			if inova_AF > 0.0:
				if inova_AF < 0.01:
				    inova_001 += 1
				elif inova_AF < 0.05:
				    inova_005 += 1
				else:
				    inova_005_plus +=1





	AF_file.close()
	filter_file.close()
	print "Total lines"
	print str(counter)
	print "New Filtered out lines"
	print str(counter_unique_no_pass)

	
	#counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/Repeats_filters_counter.txt"
	#counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/Homopolymer_filters_counter.txt"
	#counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/Microsatelite_filters_counter.txt"
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/SegDup_filters_counter.txt"
	#counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/RepeatMask_filters_counter.txt"

	counterf = open(counter_file, 'a')
	
	AF_counter = chrom + "\t" + str(counter_unique_no_pass) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
				"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	print "DONE!"




if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1], sys.argv[2])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)
