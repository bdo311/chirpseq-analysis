# getTotalCount.py
# 3/1/14
# Gets total count for all ChIP-seq reads

import csv, sys, fileinput
csv.register_dialect("textdialect", delimiter='\t')

if len(sys.argv) > 1:
	fn = sys.argv[1]

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	total = 0
	counter = 0
	for row in reader:
		#if 'chr' not in row[0]: continue
		#if counter % 100000 == 0: print counter
		total = total + (int(row[2]) - int(row[1])) * float(row[3])
		#print int(row[2])-int(row[1]),float(row[3])
		counter += 1

	print fn, 'combined', total
else:
	total = 0
	counter = 0
	for line in fileinput.input():
		row = line.split()
		total = total + (int(row[2]) - int(row[1])) * float(row[3])
		#print int(row[2])-int(row[1]),float(row[3])
		counter += 1

	print 'combined', total
		