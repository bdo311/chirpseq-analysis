# fixfastq.py
# runs through a fastq file and throws out all reads that don't make sense
# Usage: nohup cat ddx5_ActD_R1_test.fastq | parallel -j 8 --pipe -L3000 python fixfastq.py > ddx5_ActD_R1_test_ofn.fastq &


import fileinput

def writeLines(lines):
	if len(lines[1]) != len(lines[3]): return
	for l in lines:
		print l.rstrip()
		
lines = []
counter = 0
goodFlag = False
for line in fileinput.input():
	counter += 1
	if counter % 4 == 1:
		if line[0] == '@':
			if lines != [] and goodFlag: writeLines(lines)
			lines = [line]
			goodFlag = True			
		else:
			lines = []
			goodFlag = False
			counter -= 1
	else:
		lines.append(line)
if len(lines) == 4: writeLines(lines)
			
		