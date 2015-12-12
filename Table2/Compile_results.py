"""
Compile annotation files


"""
import os, sys, gzip, datetime
import numpy as np
import operator


def main():

	input_file_counts = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.final")
	input_file_novel = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.thousandG")

	outfile_file_cunts = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.compiled", "w")
	outfile_file_novel = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.thousandg.compiled", "w")
	
	count_dict = {}
	novel_dict = {}
	for line in input_file_counts:
		tp_line = line.strip().split("\t")
		dict_key = tp_line[1]
		array2 = [int(i) for i in tp_line[2:]]
		print "array 2 " + str(array2) 
		#see if this entry was already found
		try:
			array1 = count_dict[dict_key]
			print "array 1 " + str(array1)
			#add new data to the dictionary
			sum_array = [x + y for x, y in zip(array1, array2)]
			print "sum array " + str(sum_array)
			count_dict[dict_key] = sum_array
		except:
			count_dict[dict_key] = array2

	input_file_counts.close()

	for i in count_dict:
		counts = count_dict[i]
		outfile_file_cunts.write(i+ "\t" + str(counts[0])+"\t"+str(counts[1])+"\t"+str(counts[2])+"\n")
	outfile_file_cunts.close()

	#Compile the novel data
	for line in input_file_novel:
		tp_line = line.strip().split("\t")
		dict_key = tp_line[1]
		array2 = [float(i) for i in tp_line[2:]]
		#see if this entry was already found
		try:
			array1 = novel_dict[dict_key]
			sum_vector = vector1+vector2
			print "array 1 " + str(array1)
			#add new data to the dictionary
			sum_array = [x + y for x, y in zip(array1, array2)]
			print "sum array " + str(sum_array)

			#add new data to the dictionary
			novel_dict[dict_key] = sum_array
		except:
			novel_dict[dict_key] = array2

	input_file_novel.close()

	for i in novel_dict:
		counts = novel_dict[i]
		new_counts = []
		#Since we added the percentage by chromosome, we need to split by 22
		for c in counts:
			n_c = float(c)/22.0
			new_counts.append(n_c)

		outfile_file_novel.write(i+ "\t" + str(new_counts[0])+"\t"+str(new_counts[1])+"\t"+str(new_counts[2])+"\n")
	outfile_file_novel.close()




if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)