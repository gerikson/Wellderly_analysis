"""
Parse wellderly & inova from GenomeComb files

August 4th 2015

"""


import os, sys, gzip, time
from time import gmtime, strftime
import multiprocessing as mp
import subprocess
import argparse


WRITE_BUFFER_SIZE = 100
original_head_array = []
interesting_column_names = ['ploidy', 'varFilter', 'alternativeCalls', 'sequenced',
        'locus', 'totalScore', 'refcons']


def load_chrom(chrom):

    chrom_filename="/gpfs/group/stsi/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/"+chrom+".fa"
    fasta_sequences = SeqIO.parse(open(chrom_filename),'fasta')
    for fasta in fasta_sequences:
        sequence = fasta.seq.tostring()
    return sequence

#This function adds the padding for all of the non snps
def add_padding(original_line, sequence):
    begin=int(original_line[1])
    final_line = ""
    ref=original_line[4]
    alt=original_line[5]
    form="GT:VS1:VS2:VQ1:VQ2:CV:RS"
 
    if ref == '@' or ref == '':
        beg = begin - 1
        ref=sequence[beg].upper()
        if ',' in alt:
            a = alt.split(',')
            for index, val in enumerate(a):
                a[index] = sequence[beg].upper()+val
            alt = ','.join(a)
        else:
            alt=sequence[beg].upper()+alt

    #This is a insertion
    elif alt == '@' or alt == '':
        beg = begin - 1
        ref=sequence[beg].upper()+ref
        alt=sequence[beg].upper()

    elif (len(alt) > 1 or len(ref) > 1):
        beg = begin - 1
        ref=sequence[beg].upper()+ref
        if ',' in alt:
            a = alt.split(',')
            for index, val in enumerate(a):
                a[index]=sequence[beg].upper()+val
                alt = ','.join(a)
        else:
            alt=sequence[beg].upper()+alt

    else:
        print "What the heck is this? " + "\t".join(original_line)

    #add 1 from zero to 1 coordinate
    begin = begin
    final_line = original_line[0]+"\t"+str(begin)+"\t.\t" +ref+"\t"+alt+"\t.\t.\t.\t"+form

    return final_line, alt

def extract_columns(line, original_head_array):
    #print "herro from extract"
    original_line = line.split("\t")
    ref = original_line[4]

    #IF alt == '@', skip those entries
    if original_line[5] == '@':
        final_line = ""
        return final_line
    #if multiple entries, remove the '@'
    else:
        alt = original_line[5].split(',')
        if '@' in original_line[5]:
            alt.remove['@']
            original_line[5] = ",".join(alt)
            print orginal_line[5]

    final_line  = "\t".join(original_line[:6])
    
    wel_clustered = []
    inova_clustered = []
    wel_undefined = 0
    inova_undefined = 0
    #import pdb; pdb.set_trace()
    for index, col in enumerate(original_line[6:]):
        index += 6
        #extract ploidy columns
        if not (is_junk(col) or original_head_array[index][:10] == 'totalScore' \
                or original_head_array[index][:7] == 'refcons'):
            #Deal with empty columns
            if col == '':
                col='.'

            #extract zyg column
            if original_head_array[index][:3] == 'zyg':
                if col == 'r':
                    final_line=final_line+"\t0/0"
                elif col in ['c', 'm', 't', 'u', '?']:
                    wel_undefined, inova_undefined, final_line = \
                            extract_genotype(original_line, index, \
                                final_line, ref, alt, wel_undefined, \
                                inova_undefined)
                    #print final_line
                else:
                    print 'Something weird is going on:\nzyg={}\n{}'\
                            .format(col, original_head_array[index])

            elif original_head_array[index][:8] == 'varScore' or original_head_array[index][:10] == 'varQuality' \
                    or original_head_array[index][:8] == 'coverage' or original_head_array[index][:8] == 'refscore':
                    #final_line += "\t" + col
                    final_line += ":" + col

            #extract the cluster entry if any
            elif original_head_array[index][:7] == 'cluster':

                if col !="?" and col != "" and col != ".":
                    if original_head_array[index].endswith('DID'):
                        wel_clustered.append(col)
                    else:
                        inova_clustered.append(col)


    final_line = final_line.strip() + "\t" + str(wel_undefined) \
        + "\t" + str(inova_undefined) + "\t" + ",".join(wel_clustered) + "\t" + ",".join(inova_clustered) + "\n"
    
    #final_line = final_line.replace("?", ".")
    #print final_line
    return final_line
'''
#Multiprocess version of the method
def extract_columns(chunk, original_head_array):
    return_chunk = []
    for line in chunk:
        #print "herro from extract"
        original_line = line.split("\t")
        ref = original_line[4]

        #IF alt == '@', skip those entries
        if original_line[5] == '@':
            continue
        #if multiple entries, remove the '@'
        else:
            alt = original_line[5].split(',')
            if '@' in original_line[5]:
                alt.remove['@']
                original_line[5] = ",".join(alt)
                print orginal_line[5]

        final_line  = "\t".join(original_line[:6])
        
        wel_clustered = []
        inova_clustered = []
        wel_undefined = 0
        inova_undefined = 0
        #import pdb; pdb.set_trace()
        for index, col in enumerate(original_line[6:]):
            index += 6
            #extract ploidy columns
            if not (is_junk(col) or original_head_array[index][:10] == 'totalScore' \
                    or original_head_array[index][:7] == 'refcons'):
                #Deal with empty columns
                if col == '':
                    col='.'

                #extract zyg column
                if original_head_array[index][:3] == 'zyg':
                    if col == 'r':
                        final_line=final_line+"\t0/0"
                    elif col in ['c', 'm', 't', 'u', '?']:
                        wel_undefined, inova_undefined, final_line = \
                                extract_genotype(original_line, index, \
                                    final_line, ref, alt, wel_undefined, \
                                    inova_undefined)
                        #print final_line
                    else:
                        print 'Something weird is going on:\nzyg={}\n{}'\
                                .format(col, original_head_array[index])

                elif original_head_array[index][:8] == 'varScore' or original_head_array[index][:10] == 'varQuality' \
                        or original_head_array[index][:8] == 'coverage' or original_head_array[index][:8] == 'refscore':
                        #final_line += "\t" + col
                        final_line += ":" + col

                #extract the cluster entry if any
                elif original_head_array[index][:7] == 'cluster':

                    if col !="?" and col != "" and col != ".":
                        if original_head_array[index].endswith('DID'):
                            wel_clustered.append(col)
                        else:
                            inova_clustered.append(col)


        final_line = final_line.strip() + "\t" + str(wel_undefined) \
            + "\t" + str(inova_undefined) + "\t" + ",".join(wel_clustered) + "\t" + ",".join(inova_clustered)
        
        return_chunk.append(final_line)
        #final_line = final_line.replace("?", ".")
        #print final_line
    
    return return_chunk
'''

def extract_genotype(original_line, index, final_line, ref, alt, wel_undefined, inova_undefined):
	genotype = ""
	#get first alleles 

	if ref == original_line[index+1]:

		genotype = "0/"

	else:
		if original_line[index+1] == '?' or original_line[index+1] == '-':
			geno = 'N'
		else:
			geno = 'X'
		for i, j in enumerate(alt):
			if j == original_line[index+1]:
				geno = str(i+1)
		genotype = geno + "/"

	#get second allele
	if ref == original_line[index+2]:

		genotype += "0"
	else:
		if original_line[index+2] == '?' or original_line[index+2] == '-':
			geno = 'N'
		else:
			geno = 'X'
		for i, j in enumerate(alt):
			if j == original_line[index+2]:
				geno = str(i+1)

		genotype += geno

	if 'N' in genotype:
		if original_head_array[index].endswith('DID'):
			wel_undefined += 1
		else:
			inova_undefined += 1

	final_line += "\t" + genotype
	return wel_undefined, inova_undefined, final_line


def is_junk(col):
    """ Determines if the given column name is junk. """
    for name in interesting_column_names:
        if name in col:
            return True
    return False

def add_header(o):
    o.write("##fileformat=VCFv4.1\n")
    o.write("##fileDate=20150706\n")
    o.write("##source=myParsedGenomeCombFile\n")
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write('##FORMAT=<ID=VS1,Number=1,Type=String,Description="varScoreEAF1">\n')
    o.write('##FORMAT=<ID=VS2,Number=1,Type=String,Description="varScoreEAF2">\n')
    o.write('##FORMAT=<ID=VQ1,Number=1,Type=String,Description="varQuality1">\n')
    o.write('##FORMAT=<ID=VQ2,Number=1,Type=String,Description="varQuality2">\n')
    o.write('##FORMAT=<ID=CV,Number=1,Type=String,Description="coverage">\n')
    o.write('##FORMAT=<ID=RS,Number=1,Type=String,Description="refscore">\n')
    return o

def main(chrom):
    c=chrom[3:]
    #input_filename="/projects/wellderly/GenomeComb/inova_results/mcompar."+c+".txt.gz"
    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_results/mcompar."+c+".txt.gz"
    #output_filename="/projects/wellderly/GenomeComb/inova_parsed/mcompar.parced."+chrom+".txt.gz"
    output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/parsed/mcompar.parsed."+chrom+".txt.gz"
    #import pdb; pdb.set_trace()
    wellderly_total_individual = 0
    inova_total_individual = 0

    write_buffer = ""
    total_sample = 0
    block = ""
    global original_head_array
    original_head_array = []

    # Calculations
    print 'Calculating...'

    with gzip.open(input_filename) as f, gzip.open(output_filename, 'w') as o:
        #o = add_header(o)
        #sequence=load_chrom(chrom)
        #count = mp.cpu_count()
        #print str(count)
        #pool = mp.Pool(processes=(count - 3))
        #chunk = []
        counter = 0
        print strftime("%H:%M:%S", gmtime())
        for line in f:
            counter = counter + 1

            if counter == 1:
                original_head_array = line.split()
                header = "\t".join(original_head_array[:6])

                for index, col in enumerate(original_head_array[6:]):
                    index += 6
                    #extract ploidy columns and other junk columns
                    if not is_junk(col):
                        if original_head_array[index][:3] == 'zyg':
                            split_col = col.split('-')
                            header += "\tgenotype" + split_col[1] + "-" + split_col[2]
                            total_sample += 1
                            if col.endswith('DID'):
                                wellderly_total_individual += 1
                            else:
                                inova_total_individual += 1
                        elif original_head_array[index][:8] == 'varScore' or original_head_array[index][:10] == 'varQuality'\
                                or original_head_array[index][:8] == 'coverage' or original_head_array[index][:8] == 'refscore':
                            #header += "\t" + col
                            header += ":" + col

                #header += "\tMissing/Uncertain_Well\tMissing/Uncertain_Inova\tClustered_Wellderly\tClustered_Inova\tRepeatRegion\tMicrosatelite\tSegmentalDuplication\tAdjacentHomopolymer\n"
                header += "\tMissing/Uncertain_Well\tMissing/Uncertain_Inova\tClustered_Wellderly\tClustered_Inova\n"
                
                o.write(header)
                o.flush()
                print "Total Wellderly: {}\nTotal Inova: {}\nTotal Samples: {}"\
                    .format(wellderly_total_individual, inova_total_individual, total_sample)
                #print header
            
            # Parse all the non-header lines.
            else:
                
                
                final_line = extract_columns(line, original_head_array)

                block = block + final_line

                if counter%100 == 0:
                    #print strftime("%H:%M:%S", gmtime())
                    o.write(block)
                    o.flush()
                    block = ""
                
                if counter%10000 == 0:
                    print strftime("%H:%M:%S", gmtime())
                
                '''
                #multiprocess version
                chunk.append((line, original_head_array))

                if counter%100 == 0:
                    results = []
                    r = pool.map_async(extract_columns, chunk, callback=results.extend)
                    r.wait()
                    chunk = []
                    print strftime("%H:%M:%S", gmtime())
                    for res in results:
                        print res
                        o.write(res + "\n")
                        o.flush()
                    results = []    
                '''
    '''            
    results = []        
    r = pool.map_async(extract_columns, chunk, callback=results.extend)
    r.wait()
    chunk = []
    
    for res in results:
        o.write(res + "\n")
        o.flush()

    pool.close()
    pool.join()
    '''

if __name__ == '__main__':

    print "Python Version: " + sys.version

    #parser = argparse.ArgumentParser()
    #parser.add_argument('-f','--input', type=str,help='Specifies the input file, /path/to/GenomeCombData/mcompar.chrom.txt.gz')
    #parser.add_argument('-o','--output', type=str,help='Specifies the input file, /path/to/GenomeCombData/mcompar.chrom.txt.gz')
    #args = vars(parser.parse_args())
    
    start = time.time()

    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    print 'This took {0} seconds'.format(end - start)

