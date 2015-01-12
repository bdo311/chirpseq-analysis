# rescaleRepeatBedgraph.py
# 1/9/15
# reads in flagstat for repeat, and scales a normalized bedGraph

import csv, sys, math
csv.register_dialect("textdialect", delimiter='\t')

def main():
	if len(sys.argv) != 4: 
		print "Need two arguments."
		print "Usage: python rescaleRepeatBedgraph.py <flagstat> <bedgraph> <output>"
		
	fs = sys.argv[1]
	bg = sys.argv[2]
	obg = sys.argv[3]
	
	# read in flagstat	
	ifile = open(fs, 'r')
	total = ifile.readline().split(' ')[0]
	ifile.readline()
	mapped = ifile.readline().split(' ')[0]
	ifile.close()
	ratio = float(mapped)/float(total)

	# read in bedgraph
	sum = 0
	with open(bg, 'r') as ifile:	
		reader = csv.reader(ifile, 'textdialect')
		for row in reader:
			if '_' in row[0] and 'chr' in row[0]: continue
			sum += math.fabs((int(row[2]) - int(row[1])) * float(row[3]))

	# multiply ratio by 500M
	multFactor = int(5e8)/float(sum)*ratio
	with open(bg, 'r') as ifile, open(obg, 'w') as ofile:	
		reader = csv.reader(ifile, 'textdialect')
		writer = csv.writer(ofile, 'textdialect')
		for row in reader:
			if '_' in row[0] and 'chr' in row[0]: continue
			writer.writerow([row[0], row[1], row[2], multFactor*float(row[3])])
			


if __name__ == '__main__':
	main()