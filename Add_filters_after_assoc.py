"""
Add filters to the association file

Ex association file:
 CHR  SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  22    .   16050007    G        0 0.0003318    A       0.3989       0.5277            0 
  22    .   16050017    T        0 0.000332    G       0.3992       0.5275            0 
  22    .   16050021    C        0 0.0003318    T       0.3989       0.5277            0 


Adding filters at the end:

NA - means variants found in region, filter out (used NA as it's the easiest to filter out with R is.na)
G - good variant, keep
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
	#input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/association/"+chrom+".assoc"
	repeats_filename="/gpfs/home/gerikson/wellderly/filter_data/SimpleRepeats/byChrom/simpleRepeats."+chrom+".txt"
	micrositelite_filename="/gpfs/home/gerikson/wellderly/filter_data/Micrositelite/byChrom/micrositelite."+chrom+".txt"
	homopolymers_filename="/gpfs/home/gerikson/wellderly/filter_data/Homopolymer/byChrom/homopolymer."+chrom+".txt"
	segmentalDup_filename="/gpfs/home/gerikson/wellderly/filter_data/SegmentalDup/byChrom/SegmentalDup."+chrom+".txt"
	RepeatMasker_filename="/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom_ALL/RepeatMasker."+chrom+".txt"
	#output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/association_with_filters"+chrom+".assoc"

	#For testing
	#input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/bed_files/chr22.0.01AF.assoc"
	#output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/bed_files/chr22.0.01AF.withFilters.assoc"
	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/plink/chr22_VQHIGH0.95_excludeMissingGeno0.1/chr22.0.95.snp.noMissingGeno.maf.assoc"
	output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/plink/chr22_VQHIGH0.95_excludeMissingGeno0.1/chr22.withFilters0.95snp.maf.assoc"
	repeats_begin = []
	micrositelites_begin = []
	homopolymers_begin = []
	segmentalDup_begin = []
	repeatMask_begin = []

	repeats_end = []
	micrositelites_end = []
	homopolymers_end = []
	segmentalDup_end = []
	repeatMask_end = []

	# Prep
	print 'Prepping repeats...'
	with open(repeats_filename) as f:
	    for line in f:
	        items = line.split('\t')
	        repeats_begin.append(int(items[2]))
	        repeats_end.append(int(items[3]))

	
	# Micrositelites
	print 'Prepping micrositelites...'
	with open(micrositelite_filename) as f:
	    for line in f:
			items = line.split('\t')
			micrositelites_begin.append(int(items[2]))
			micrositelites_end.append(int(items[3]))
	#Homopolymer
	print 'Prepping homopolymers...'
	with open(homopolymers_filename) as f:
	    for line in f:
	        items = line.split('\t')
	        homopolymers_begin.append(int(items[1]))
	        homopolymers_end.append(int(items[2]))

	#Segmental Duplication
	print 'Prepping SegmentalDup...'
	with open(segmentalDup_filename) as f:
	    for line in f:
	        items = line.split('\t')
	        segmentalDup_begin.append(int(items[2]))
	        segmentalDup_end.append(int(items[3]))
	
	#Segmental Duplication
	print 'Prepping RepeatMask...'
	with open(RepeatMasker_filename) as f:
	    for line in f:
	        items = line.split('\t')
	        repeatMask_begin.append(int(items[6]))
	        repeatMask_end.append(int(items[7]))

	print 'Calculating...'

	counter = 0
	repeat_counter = 0
	microsat_counter = 0
	homo_counter = 0
	segDup_counter = 0
	repeatMask_counter = 0
	block = ""

	with open(input_filename) as f:
		with open(output_filename, 'w') as o: 

			for line in f:
				counter += 1
				# Parse all the non-header lines.
				if counter == 1:
					line = line.strip().split()
					o.write("\t".join(line) +"\tRepeat\tMicrositelite\tHomopolymer\tSegmentalDup\tRepeatMasker\n")
				else:

					#print line
					line = line.strip().split()
					begin = line[2]
					#try:
					end = int(line[2]) + len(line[3])
					#except:
					#	print line
					final_line = "\t".join(line)
					#ADD repeats
					repeatOverlap = Overlap(begin, repeats_begin, end, repeats_end)
					overlapRep = repeatOverlap.main()
					if overlapRep:
						repeat_counter += 1
						final_line += "\tNA"
					else:
						final_line +="\tG"

					#ADD micrositelites
					microOverlap = Overlap(begin, micrositelites_begin, end, micrositelites_end)
					overlapMicro = microOverlap.main()
					if overlapMicro:
						microsat_counter += 1
						final_line += "\tNA"
					else:
						final_line +="\tG"

					#Add homopolymers
					homoOverlap = Overlap(begin, homopolymers_begin, end, homopolymers_end)
					overlapHomo = homoOverlap.main()
					if overlapHomo:
						homo_counter += 1
						final_line += "\tNA"
					else:
						final_line +="\tG"


					#Add SegmentalDup
					segOverlap = Overlap(begin, segmentalDup_begin, end, segmentalDup_end)
					overlapSeg = segOverlap.main()
					if overlapSeg:
						segDup_counter += 1
						final_line += "\tNA\t"
					else:
						final_line +="\tG\t"

					#Add RepeatMask
					repOverlap = Overlap(begin, repeatMask_begin, end, repeatMask_end)
					overlapRep = repOverlap.main()
					if overlapRep:
						repeatMask_counter += 1
						final_line += "\tNA\n"
					else:
						final_line +="\tG\n"

					block = block + final_line
					if counter%100 == 0:
						o.write(block)
						block = ""
						o.flush()

					if counter%10000 == 0:
						print datetime.datetime.now().time()
						print str(counter)
						sys.stdout.flush()

			o.write(block)
	block = ""
	print "Total lines " + str(counter)
	print "Total repeats counter " + str(repeat_counter)
	print "Total microsat_counter " + str(microsat_counter)
	print "Total homo counter " + str(homo_counter)
	print "Total segDup counter " + str(segDup_counter)
	print "Total repeatMask counter " + str(repeatMask_counter)


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)