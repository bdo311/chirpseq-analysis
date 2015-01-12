# countRepeatReads.py
# 12/13/14

import csv, os, sys, multiprocessing, glob, subprocess
csv.register_dialect("textdialect", delimiter='\t')

def getStats(fn, d, regionToLocus):
	nums = []
	cmd = "samtools view -F 4 " + fn + " | wc -l"
	result = subprocess.check_output(cmd, shell=True).rstrip()
	nums.append(result) #number of mapped reads
	
	for region in regionToLocus.keys():
		locus = regionToLocus[region]
		cmd = "samtools view " + fn + " " + locus + " | wc -l"
		result = subprocess.check_output(cmd, shell=True).rstrip()
		print cmd, result
		nums.append(result)		
	d[fn[:-4]] = nums
	return d
	
def main():
	if len(sys.argv)<4: 
		print "Usage: python countMappedReads.py <config file> <position file> <output file>"
		exit()
		
	conf = sys.argv[1]
	pos = sys.argv[2]
	ofn = sys.argv[3]

	# read list of files that need to be read
	files = []
	ifile = open(conf, 'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:
		if len(row) > 0: files.append(row[0])
	ifile.close()
	
	# read position file
	ifile = open(pos, 'r')
	reader = csv.reader(ifile, 'textdialect')
	regionToLocus = {}
	for row in reader:
		regionToLocus[row[0]] = 'Repeat_regions:' + row[2] + '-' + row[3]
	ifile.close()
	
	# get output file
	d = {}
	for fn in files:
		getStats(fn, d, regionToLocus)
	
	# write into output file
	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	header = ['', 'mapped']
	header.extend(regionToLocus.keys())
	writer.writerow(header)
	for fn in d:
		outputRow = [fn]
		outputRow.extend(d[fn])
		writer.writerow(outputRow)
		
	ofile.close()
	
if __name__ == "__main__":
	main()
