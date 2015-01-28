# makeBedgraphSmaller.py
# 8/5/14
# compresses a bedgraph by making the step size >= a certain threshold

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')

def main():
	ifn, ofn = sys.argv[1], sys.argv[2]
	thresh = int(sys.argv[3])

	ifile = open(ifn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')

	# variables that carry forward
	currStart = 0
	currLength = 0
	currTotal = 0

	# fencepost first row
	row = reader.next()
	newTranscript = True
	regionStart, regionEnd = int(row[1]), int(row[2])
	regionLength = regionEnd - regionStart
	regionTotal = float(row[3]) * regionLength
	if regionLength >= thresh: writer.writerow(row)
	else: 
		currStart = regionStart
		currLength = regionLength
		currTotal = regionTotal
		newTranscript = False

	# all other rows
	for row in reader:
		regionStart, regionEnd = int(row[1]), int(row[2])
		regionLength = regionEnd - regionStart
		regionTotal = float(row[3]) * regionLength

		if newTranscript: 
			if regionLength >= thresh: writer.writerow(row)
			else: 
				currStart = regionStart
				currLength = regionLength
				currTotal = regionTotal
				newTranscript = False
		else: 
			currLength = regionEnd - currStart
			currTotal = currTotal + regionTotal
			if currLength >= thresh: 
				writer.writerow([row[0], currStart, regionEnd, currTotal/currLength])
				newTranscript = True

	ifile.close()
	ofile.close()

if __name__ == '__main__':
	main()
