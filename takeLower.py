# takeLower.py
# 6/16/14
# takes the lower of even and odd in order to merge the bedgraph

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')

ifn = sys.argv[1]
ofn = sys.argv[2]

ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

reader.next()
for row in reader:
	num = min(float(row[3]), float(row[4]))
	if num == 0: continue
	outputRow = row[:3]
	outputRow.append(num)
	writer.writerow(outputRow)

ifile.close()
ofile.close()
