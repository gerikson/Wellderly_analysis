"""
Remove the variants that fall into any of these regions:
Reapeats
Homopolymers
Micrositelites
segmentalDup
RepeatMasker

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
 


def main(chrom):

	repeats_filename="/gpfs/home/gerikson/wellderly/filter_data/SimpleRepeats/byChrom/simpleRepeats."+chrom+".txt"
	micrositelite_filename="/gpfs/home/gerikson/wellderly/filter_data/Micrositelite/byChrom/micrositelite."+chrom+".txt"
	homopolymers_filename="/gpfs/home/gerikson/wellderly/filter_data/Homopolymer/byChrom/homopolymer."+chrom+".txt"
	segmentalDup_filename="/gpfs/home/gerikson/wellderly/filter_data/SegmentalDup/byChrom/SegmentalDup."+chrom+".txt"
	RepeatMasker_filename="/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom/RepeatMasker."+chrom+".txt"

	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered/wellderly_inova.VQHIGH.0.95white.nocluster."+str(chrom)+".vcf.gz"
	output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc/wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc."+str(chrom)+".vcf.gz"
	counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc/counter_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.txt"

	#filter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered/Filter_file_repeat_homopoly_etc."+str(chrom)+".txt"
	#filt = open(filter_file, "w")
	#filt.write("#CHROM\tPOS\tID\tREF\tALT\tRepeat\tMicrositelite\tHomopolymer\tSegmentalDul\tReapeatMask\n")

	c = open(counter_file, "a")
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
	good_lines = 0

	repeat_counter = 0
	microsat_counter = 0
	homo_counter = 0
	segDup_counter = 0
	repeatMask_counter = 0
	block = ""
	f = gzip.open(input_filename)
	o = gzip.open(output_filename, 'w')
	#filter_block = ""

	for line in f:

		if line[:1] == "#":
			print "header"
			o.write(line)

		else:

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


			tp_line = line.strip().split()
			begin = tp_line[1]

			end = int(int(tp_line[1])) + int(len(tp_line[3]))


			#filter_block = filter_block + "\t".join(tp_line[:6])

			#ADD repeats
			repeatOverlap = Overlap(begin, repeats_begin, end, repeats_end)
			overlapRep = repeatOverlap.main()
			if overlapRep:
				repeat_counter += 1
				#filter_block = filter_block + "\tYES\t\t\t\t\n"
				continue
			else:
				block = block + line
				good_lines += 1
				#filter_block = filter_block + "\t"
			

			#ADD micrositelites
			microOverlap = Overlap(begin, micrositelites_begin, end, micrositelites_end)
			overlapMicro = microOverlap.main()
			if overlapMicro:
				microsat_counter += 1
				#filter_block = filter_block + "\tYES\t\t\t\n"
				continue
			
			else:
				block = block + line
				good_lines += 1
				#filter_block = filter_block + "\t"
			
			#Add homopolymers
			homoOverlap = Overlap(begin, homopolymers_begin, end, homopolymers_end)
			overlapHomo = homoOverlap.main()
			if overlapHomo:
				homo_counter += 1
				#filter_block = filter_block + "\tYES\t\t\n"
				continue
			else:
				block = block + line
				good_lines += 1
				#filter_block = filter_block + "\t"


			#Add SegmentalDup
			segOverlap = Overlap(begin, segmentalDup_begin, end, segmentalDup_end)
			overlapSeg = segOverlap.main()
			if overlapSeg:
				segDup_counter += 1
				#filter_block = filter_block + "\tYES\t\n"
				continue
			else:
				#filter_block = filter_block + "\t"
				block = block + line
				good_lines += 1

			#Add RepeatMask
			repOverlap = Overlap(begin, repeatMask_begin, end, repeatMask_end)
			overlapRep = repOverlap.main()
			if overlapRep:
				repeatMask_counter += 1
				#filter_block = filter_block + "\tYES\n"
				continue
			else:
				#filter_block = filter_block + "\t\n"
				block = block + line
				good_lines += 1

			

	#filt.write(filter_block)
	o.write(block)
	print "Total repeats counter " + str(repeat_counter)
	print "Total microsat_counter " + str(microsat_counter)
	print "Total homo counter " + str(homo_counter)
	print "Total segDup counter " + str(segDup_counter)
	print "Total repeatMask counter " + str(repeatMask_counter)
	
	print "Total lines " + str(counter)
	print "Total good lines " + str(good_lines)

	c.write(chrom + "\t" + str(counter) + "\t" + str(good_lines)+ "\n")
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