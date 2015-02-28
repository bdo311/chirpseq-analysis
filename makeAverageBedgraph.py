# makeAverageBedgraph.py
# 2/27/15
# averages multiple bedgraphs made from unionbedg

import csv, sys
import numpy as np
csv.register_dialect("textdialect", delimiter='\t')

ifn = sys.argv[1]
ofn = sys.argv[2]

ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	num = np.mean([float(x) for x in row[3:]])
	if num == 0: continue
	outputRow = row[:3]
	outputRow.append(num)
	writer.writerow(outputRow)

ifile.close()
ofile.close()
