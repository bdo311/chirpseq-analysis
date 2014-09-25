# normalizeBedgraph.py
# 6/16/14
# normalizes one bedgraph at a time

import csv, sys, math
csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("spacedialect", delimiter=' ')

ifn = sys.argv[1]
sizes = sys.argv[2]
ofn = sys.argv[3]

# get biggest allowed sequence sizes
ifile = open(sizes, 'r')
reader = csv.reader(ifile, 'textdialect')

chrToLimit = {}
for row in reader:
	chrToLimit[row[0]] = int(row[1])

ifile.close()

# normalize by scaling everything to a total read density of 500 million
ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')

sum = 0
for row in reader:
	if '_' in row[0]: continue
	sum += math.fabs((int(row[2]) - int(row[1])) * float(row[3]))

multFactor = 500000000/sum

ifile.seek(0)

ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

pastRow = ['0', 0, 0, 0]
for row in reader:
	if '_' in row[0]: continue
	start, stop, val = int(row[1]), int(row[2]), float(row[-1])
	if val == 0: continue
	if start > chrToLimit[row[0]] or stop > chrToLimit[row[0]]: continue
	val = val * multFactor
	row[-1] = val
	
	if start == int(pastRow[2]) and val == float(pastRow[-1]): 
		pastRow[2] = stop
	else:
		if pastRow[0] != '0': writer.writerow(pastRow)
		pastRow = row
	
	
if pastRow[0] != '0': writer.writerow(pastRow)
ifile.close()
ofile.close()

