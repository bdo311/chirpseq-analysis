# runRetroviral.py
# 3/27/2016
# Makes histograms and count tables for retroviral data; requires either a
# trimmed FASTQ file or a STAR-generated BAM file as input

import sys
import os
import argparse
import pandas as pd
import numpy as np
import collections
import csv
csv.register_dialect("textdialect", delimiter='\t')

### Parsing arguments ###
parser = argparse.ArgumentParser(description="runRetroviral.py", epilog="Example: python runRetroviral.py ")
parser.add_argument('-i', metavar='INPUT', help="Input fastq or bam file", required=True)
group = parser.add_mutually_exclusive_group()
group.add_argument('--fastq', action='store_true', help="required if your input is a trimmed FASTQ or FASTQ.GZ")
group.add_argument('--bam', action='store_true', help="required if your input is a trimmed BAM")
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('--human', action='store_true', help="required if your file is from human")
group2.add_argument('--mouse', action='store_true', help="required if your file is from mouse")
parser.add_argument('--compressed', action='store_true', help="required if your input is *.gz")
parser.add_argument('--iclip', action='store_true', help="specify if file is from CLIP data")
parser.add_argument('-p', metavar='PREFIX', help="Prefix of output files", required=True)
parser.add_argument('-sr', metavar='STAR_RATIO', type=float, help="Maximum mismatches per base allowed for STAR genome mapping (corresponds to outFilterMismatchNoverLmax). Default is 0.08 (2 mismatches per 25 mapped bases).", default=0.08)

# organism
args = parser.parse_args()
if not (args.human or args.mouse):
	print "Error: must include --human or --mouse. Exiting."
	exit()
if args.human: org = 'human' 
else: org = 'mouse'

# filetype
if not (args.fastq or args.bam):
	print "Error: must include --fastq or --bam. Exiting."
	exit()
if args.fastq: filetype = 'fastq' 
else: filetype = 'bam'

# other arguments
compressed = True if args.compressed else False
iclip = True if args.iclip else False
ifile = args.i
prefix = args.p
star_ratio = args.sr

home = "/home/raflynn/CLIP/"
if org == 'human':
	repeat_index = home + '/docs/hg19/repeat/rep_spaced' # bt2 index for repeat RNA.
	retro_index = home + '/docs/hg19/retroviral/'
	retro_annotation = retro_index + 'Hs_retroviralIndex_spaced_positions.txt'
elif org == 'mouse':
	repeat_index = home + '/docs/mm9/repeat/rep_spaced' # bt2 index for repeat RNA.
	retro_index = home + '/docs/mm9/retroviral/'
	retro_annotation = retro_index + 'Mm_retroviralIndex_spaced_positions.txt'

def make_bam(ifile):

	# map to repeat
	print "Mapping fastq to repeat"
	rep_unmapped = prefix + "_repeat_unmapped.fastq"
	if compressed or ifile[-3:] == '.gz':
		cmd = "zcat {} | bowtie2 -p 8 -x {} - --un {} > /dev/null".format(ifile, repeat_index, rep_unmapped)
	else:
		cmd = "cat {} | bowtie2 -p 8 -x {} - --un {} > /dev/null".format(ifile, repeat_index, rep_unmapped)
	os.system(cmd)		
	
	# map to retro
	print "Mapping fastq to retroviral"
	retro_prefix = prefix + "_retro_"
	retro_mapped = prefix + "_retro_mapped.sam"
	cmd = "STAR --genomeDir {} --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn {} --outFileNamePrefix {} --alignEndsType EndToEnd --outFilterMismatchNoverLmax {}".format(retro_index, rep_unmapped, retro_prefix, star_ratio)
	os.system(cmd)
	os.system("mv {} {}".format(retro_prefix + "Aligned.out.sam", retro_mapped))
	
	# make bam file
	print "Making BAM file"
	bamfile_sort = prefix + "_retro_mapped.bam"
	cmd = "cat {} | samtools view -Su - -o - | samtools sort - {}".format(retro_mapped, bamfile_sort)
	os.system(cmd)
	
	return bamfile_sort

def make_hist(bed):
	data = pd.read_csv(bed, sep='\t', names=['name', 'start', 'end', 'read', 'qual', 'strand'])	
	nums = []
	for i in data.iterrows():
		nums.extend(range(i[1].start, i[1].end))
	h = np.histogram(nums, bins=np.max(nums)-np.min(nums))[0]
	hist = list(np.zeros(np.min(nums)))
	hist.extend(list(h))
	return np.array(hist)	
	
def print_list_of_hist(hist, retro):
	ofn = prefix + "_retro_hist.txt"
	with open(ofn, 'w') as ofile:
		writer = csv.writer(ofile, 'textdialect')

		numBins = 100
		for region in retro.iterrows():
			row = [region[1][0] + '__' + str(region[1].len)]
			
			nums = hist[region[1].start:region[1].end]
			intended_len = region[1].end - region[1].start
			if len(nums) != intended_len:
				nums = list(nums)
				nums.extend([0]*(intended_len - len(nums)))
				nums = np.array(nums)
			
			a = map(lambda x: float(x)*len(nums)/numBins, range(numBins)) # convert my desired scale to the current scale
			new_nums = np.interp(a, range(len(nums)), nums)
			row.extend(list(new_nums))
			
			writer.writerow(row)
			
def print_list_of_counts(bed, retro):
	data = pd.read_csv(bed, sep='\t', names=['name', 'start', 'end', 'read', 'qual', 'strand'])
	
	ofn = prefix + "_retro_cts.txt"
	fam_to_count = collections.defaultdict(lambda: 0)
	with open(ofn, 'w') as ofile:
		writer = csv.writer(ofile, 'textdialect')

		for region in retro.iterrows():
			name = region[1][0] + '__' + str(region[1].len)
			fam = name.split('__')[1]
			ct = data[data.start >= region[1].start][data.start < region[1].end].shape[0]
			fam_to_count[fam] += ct
			row = [name, region[1].len, ct]
			writer.writerow(row)
			
	ofn = prefix + "_retro_family_cts.txt"	
	with open(ofn, 'w') as ofile:
		writer = csv.writer(ofile, 'textdialect')
		for f in fam_to_count:
			writer.writerow([f, fam_to_count[f]])
		
def main():

	# read in annotation
	retro = pd.read_csv(retro_annotation, sep='\t', names=['name', 'len', 'start', 'end'])

	# 0a. if FASTQ, then map to repeat, then map to retro, then make multi BAM
	if filetype == 'fastq':
		bam = make_bam(ifile)
	else:
		bam = ifile
		
	# 1. with multi BAM make a histogram
	# print "Making histograms for each retro element"
	# bed = prefix + "_retro_mapped.bed"
	# cmd = "bamToBed -i {} > {}".format(bam, bed)
	# os.system(cmd)
	
	# total_hist = make_hist(bed)
	# print_list_of_hist(total_hist, retro)
	
	# 2. make a uni BAM
	unibam = prefix + "_retro_mapped_uni.bam"
	cmd = "samtools view -b -F 0x100 {} -o {}".format(bam, unibam)
	os.system(cmd)
	
	# 3. with uni BAM make a file of counts (by family and by subfamily)
	print "Counting reads for each retro element"

	unibed = prefix + "_retro_mapped_uni.bed"
	cmd = "bamToBed -i {} > {}".format(unibam, unibed)
	os.system(cmd)
	
	print_list_of_counts(unibed, retro)
	

if __name__ == "__main__":
	main()