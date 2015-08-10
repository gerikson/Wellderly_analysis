"""
Parse wellderly & inova from GenomeComb files

If it's deletion or insertion add padding
Else add +1 to the beggin position
Exclude the variants and the alleles that are not present in the wellderly dataset
Exclude extra columns
Format files into a vcf format

Needs module load python/2.6.5 for SeqIO

August 4th 2015

* Not properly sorted (most likely due to anchor base insertion on indels).
Will do later
* Malformed files (double tab after the FORMAT column).
Fixed
* Malformed files (invalid "@" alleles declared in ALT field and referenced by samples, 225 cases on chr1).
Fixed
* Declared ALT alleles not used (414313 cases on chr1).
Fixed
* Contig and reference lines missing (recommended by samtools validate).
Fixed
* Use of "?" to denote missing values (should really be ".").
Fixed
"""


import os, sys, gzip, time
from time import gmtime, strftime
import multiprocessing as mp
import subprocess
import argparse
from Bio import SeqIO
import datetime
#from Bio.SeqIO import FastaIO

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
def add_padding(original_line, sequence, alt_array):
    begin=int(original_line[1])
    final_line = ""
    ref=original_line[4]
    #alt=original_line[5]
    alt = ",".join(alt_array)
    form="GT:VS1:VS2:VQ1:VQ2:CV:RS"
 
    if ref == '':
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
    elif alt == '':
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


    final_line = original_line[0]+"\t"+str(begin)+"\t.\t" +ref+"\t"+alt+"\t.\t.\t.\t"+form

    return final_line

def extract_columns(line, original_head_array, sequence):
    #print "herro from extract"
    original_line = line.split("\t")
    ref = original_line[4]
    final_line = ""
    form="GT:VS1:VS2:VQ1:VQ2:CV:RS"
    #IF alt == '@', skip those entries
    if original_line[5] == '@':
        final_line = ""
        return final_line
    #if multiple entries, remove the '@'
    else:
        alt = original_line[5].split(',')
        if '@' in original_line[5]:
            alt.remove('@')
            original_line[5] = ",".join(alt)
            if original_line[5] == "":
                return ""
            #print original_line[5]

    #fix the alt allele (only when there are alternate allele)
    if original_line[5] != "":
        alt = fix_alt_allele(original_line, original_head_array)
        if len(alt) == 0:
            #Nothing to return, None of the alternate alleles are present in the wellderly
            #print "no allele found, skip this shit" 
            #print original_line[:100]
            return ""
        else:
            #original_line[5] = ",".join(alt)

            #If this is a insertion add the padding
            if ref == "":
                #print "This is a insertion "
                #print original_line[:6]
                begin=int(original_line[1])
                final_line = add_padding(original_line, sequence, alt)
                #format genotypes
                original_line[5] = ",".join(alt)
                final_geno = format_genotypes(original_line, original_head_array)
                final_line = final_line.strip() + "\t" + final_geno.strip() + "\n"
                #print "Insertion found!"
                #print final_geno
            # this is snp or sub
            else:
                begin=int(original_line[1]) + 1
                original_line[5] = ",".join(alt)
                final_line = original_line[0]+"\t"+str(begin)+"\t.\t" +ref+"\t"+','.join(alt)+"\t.\t.\t.\t"+form
                final_geno = format_genotypes(original_line, original_head_array)
                final_line = final_line.strip() + "\t" + final_geno.strip() + "\n"
                #print "Normal Variant found!"
                #print final_geno

    else:
        # If the alternate allele was found
        if find_alt_for_deletions(original_line, original_head_array):
            begin=int(original_line[1])
            final_line = add_padding(original_line, sequence, alt)
            #format genotypes
            final_geno = format_genotypes(original_line, original_head_array)
            final_line = final_line.strip() + "\t" + final_geno.strip() + "\n"
            #print "Deletion found!"
            #print final_geno
        else:
            return ""
    #print final_line
    return final_line


def find_alt_for_deletions(original_line, original_head_array):
    for index, col in enumerate(original_line[6:9075]):
        index += 6

        #extract zyg column
        if original_head_array[index][:3] == 'zyg':

            if col in ['c', 'm', 't', 'u', '?']:

                if original_line[index+1] == "" or original_line[index+2] == "":
                    return True

    #If True wasn't return that means the deletion wasn't found
    return False

#Extract the alternate alleles that are actually present in the dataset
def fix_alt_allele(original_line, original_head_array):

    #print "Begin alt set "
    #print original_line[5]

    alt = original_line[5].split(',')
    alt_set = set(alt)

    real_allele_set = set()
    for index, col in enumerate(original_line[6:9075]):
        index += 6

        #extract zyg column
        if original_head_array[index][:3] == 'zyg':

            if col in ['c', 'm', 't', 'u', '?']:
                real_allele_set.add(original_line[index+1])
                real_allele_set.add(original_line[index+2])

    #Extract the intersection of the alt set and real alt set
    final_alt = alt_set & real_allele_set

    alt = list(final_alt)
    
    #Sort the list by the lenght of alternate allele (not alphabetical)
    alt.sort(key = len)
    #print "End alt set " 
    #print alt
    return alt

#Find the genotype information, get rid of the extra data
def format_genotypes(original_line, original_head_array):
    final_line = ""
    alt = original_line[5].split(',')
    for index, col in enumerate(original_line[6:9075]):
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
                    final_line = extract_genotype(original_line, index, \
                                final_line, original_line[4], alt)
                    #print final_line
                else:
                    print 'Something weird is going on:\nzyg={}\n{}'\
                            .format(col, original_head_array[index])

            elif original_head_array[index][:8] == 'varScore' or original_head_array[index][:10] == 'varQuality' \
                    or original_head_array[index][:8] == 'coverage' or original_head_array[index][:8] == 'refscore':
                    #final_line += "\t" + col
                    final_line += ":" + col

    final_line = final_line.replace('?','.')
    return final_line

#Format the genotypes: 0, 1, . based on the vcf file format
def extract_genotype(original_line, index, final_line, ref, alt):

	genotype = ""
	#get first alleles 
	if ref == original_line[index+1]:

		genotype = "0/"

	else:
		geno = '.'
		for i, j in enumerate(alt):
			if j == original_line[index+1]:
				geno = str(i+1)
		genotype = geno + "/"

	#get second allele
	if ref == original_line[index+2]:

		genotype += "0"
	else:
		geno = '.'
		for i, j in enumerate(alt):
			if j == original_line[index+2]:
				geno = str(i+1)

		genotype += geno

	final_line += "\t" + genotype
	return final_line


def is_junk(col):
    """ Determines if the given column name is junk. """
    for name in interesting_column_names:
        if name in col:
            return True
    return False

def main(chrom):
    c=chrom[3:]
    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_results/mcompar."+c+".txt.gz"
    output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_601CG_welderly/wellderly_601CG_"+chrom+".txt.gz"
    wellderly_total_individual = 0
    inova_total_individual = 0

    write_buffer = ""
    total_sample = 0
    block = ""
    global original_head_array
    original_head_array = []

    f = gzip.open(input_filename)
    o = gzip.open(output_filename, 'w')
    # Calculations
    print 'Calculating...'
    chrom_filename="##reference=file:///gpfs/group/stsi/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/"+chrom+".fa\n"
    #with gzip.open(input_filename) as f:
    #    with gzip.open(output_filename, 'w') as o:
    o.write("##fileformat=VCFv4.1\n")
    o.write("##fileDate=20150807\n")
    o.write("##source=GenomeCombFile\n")
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write('##FORMAT=<ID=VS1,Number=1,Type=String,Description="varScoreEAF1">\n')
    o.write('##FORMAT=<ID=VS2,Number=1,Type=String,Description="varScoreEAF2">\n')
    o.write('##FORMAT=<ID=VQ1,Number=1,Type=String,Description="varQuality1">\n')
    o.write('##FORMAT=<ID=VQ2,Number=1,Type=String,Description="varQuality2">\n')
    o.write('##FORMAT=<ID=CV,Number=1,Type=String,Description="coverage">\n')
    o.write('##FORMAT=<ID=RS,Number=1,Type=String,Description="refscore">\n')
    o.write('##contig=<ID=M,length=16571,assembly=hg19>\n')
    o.write('##contig=<ID=1,length=249250621,assembly=hg19>\n')
    o.write('##contig=<ID=2,length=243199373,assembly=hg19>\n')
    o.write('##contig=<ID=3,length=198022430,assembly=hg19>\n')
    o.write('##contig=<ID=4,length=191154276,assembly=hg19>\n')
    o.write('##contig=<ID=5,length=180915260,assembly=hg19>\n')
    o.write('##contig=<ID=6,length=171115067,assembly=hg19>\n')
    o.write('##contig=<ID=7,length=159138663,assembly=hg19>\n')
    o.write('##contig=<ID=8,length=146364022,assembly=hg19>\n')
    o.write('##contig=<ID=9,length=141213431,assembly=hg19>\n')
    o.write('##contig=<ID=10,length=135534747,assembly=hg19>\n')
    o.write('##contig=<ID=11,length=135006516,assembly=hg19>\n')
    o.write('##contig=<ID=12,length=133851895,assembly=hg19>\n')
    o.write('##contig=<ID=13,length=115169878,assembly=hg19>\n')
    o.write('##contig=<ID=14,length=107349540,assembly=hg19>\n')
    o.write('##contig=<ID=15,length=102531392,assembly=hg19>\n')
    o.write('##contig=<ID=16,length=90354753,assembly=hg19>\n')
    o.write('##contig=<ID=17,length=81195210,assembly=hg19>\n')
    o.write('##contig=<ID=18,length=78077248,assembly=hg19>\n')
    o.write('##contig=<ID=19,length=59128983,assembly=hg19>\n')
    o.write('##contig=<ID=20,length=63025520,assembly=hg19>\n')
    o.write('##contig=<ID=21,length=48129895,assembly=hg19>\n')
    o.write('##contig=<ID=22,length=51304566,assembly=hg19>\n')
    o.write('##contig=<ID=X,length=155270560,assembly=hg19>\n')
    o.write('##contig=<ID=Y,length=59373566,assembly=hg19>\n')
    o.write(chrom_filename)
    counter = 0
    excluded_lines = 0
    sequence=load_chrom(chrom)
    #print strftime("%H:%M:%S", gmtime())
    for line in f:
        line = line.strip()
        counter = counter + 1

        if counter == 1:
            original_head_array = line.split("\t")
            header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

            for index, col in enumerate(original_head_array):
                #index += 6
                #extract ploidy columns and other junk columns
                if not is_junk(col) and original_head_array[index].endswith("DID"):
                    if original_head_array[index][:3] == 'zyg':
                        split_col = col.split('-')
                        header += "\t" + split_col[1] + "-" + split_col[2]
                        total_sample += 1
                        if col.endswith('DID'):
                            wellderly_total_individual += 1
                        else:
                            inova_total_individual += 1
                    elif original_head_array[index][:8] == 'varScore' or original_head_array[index][:10] == 'varQuality'\
                            or original_head_array[index][:8] == 'coverage' or original_head_array[index][:8] == 'refscore':
                        #header += "\t" + col
                        #header += ":" + col
                        continue
                elif original_head_array[index].endswith("ASM"):
                    print "Inova starts " + str(index)
                    break

            #header += "\tMissing/Uncertain_Well\tMissing/Uncertain_Inova\tClustered_Wellderly\tClustered_Inova\tRepeatRegion\tMicrosatelite\tSegmentalDuplication\tAdjacentHomopolymer\n"
            #header += "\tMissing/Uncertain_Well\tMissing/Uncertain_Inova\tClustered_Wellderly\tClustered_Inova\n"
            header = header+"\n"
            o.write(header)
            o.flush()
            #print "Total Wellderly: {}\nTotal Inova: {}\nTotal Samples: {}"\
            #    .format(wellderly_total_individual, inova_total_individual, total_sample)
            #print header
        
        # Parse all the non-header lines.
        else:
            
            
            final_line = extract_columns(line, original_head_array, sequence)

            block = block + final_line

            if final_line == "":
                excluded_lines = excluded_lines + 1
            if counter%100 == 0:
                #print strftime("%H:%M:%S", gmtime())
                o.write(block)
                o.flush()
                block = ""
            
            if counter%10000 == 0:
                print datetime.datetime.now().time()
                print str(counter)
                #print strftime("%H:%M:%S", gmtime())

    o.write(block)
    print "Total lines parsed: " + str(counter)
    #ch_count.write("chrom\ttotal_lines\texcluded_lines\n")
    chrom_count = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_601CG_welderly/chrom_count.txt"
    ch_count = open(chrom_count, 'w')
    ch = chrom +"\t"+str(counter) + "\t" + str(excluded_lines) + "\n"
    ch_count.write(ch)
    ch_count.close()
    o.flush()
    o.close()
    f.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    start = time.time()

    main(sys.argv[1])

    end = time.time()
    print "Done!"

