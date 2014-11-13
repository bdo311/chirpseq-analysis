# removeInterval.py
# 9/23/14
# zeros out the bedgraph for all intervals in the bed file; expands intervals by 2x

import csv, sys, collections
csv.register_dialect('textdialect', delimiter='\t')

def processBedFile(bed):
	ifile = open(bed, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	chrToBedInts = collections.defaultdict(lambda: [])
	for row in reader:
		start, end = int(row[1]), int(row[2])
		intLength = end - start + 1
		start = start - intLength/2
		end = end + intLength/2
		chrToBedInts[row[0]].append([start, end])
	
	ifile.close()
	return chrToBedInts

def processBedgraph(bg, chrToBedInts, ofn):
	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	ifile = open(bg, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	reader.next()
	for row in reader:
		chr = row[0]
		start, end = int(row[1]), int(row[2])
		val = float(row[3])
		
		inInterval = False
		for interval in chrToBedInts[chr]:
			if start > interval[1] or end < interval[0]: # region doesn't overlap interval at all
				continue
			elif start < interval[0]: 
				writer.writerow([chr, start, interval[0], val]) #region overlapping left side of interval
				inInterval = True
				break
			elif end > interval[1]: 
				writer.writerow([chr, interval[1], end, val])
				inInterval = True
				break
			else: #region completely inside interval
				inInterval = True
				break 
		if not inInterval: writer.writerow(row)
			
	ifile.close()
	ofile.close()

def main():
	if len(sys.argv) < 4: 
		print "Usage: python removeInterval.py <bedgraph> <bed> <ofn>"
		exit()
	bg = sys.argv[1]
	bed = sys.argv[2]
	ofn = sys.argv[3]
	
	chrToBedInts = processBedFile(bed)
	processBedgraph(bg, chrToBedInts, ofn)

if __name__ == "__main__":
	main()