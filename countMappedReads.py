# countMappedReads.py
# 12/13/14

import csv, os, sys, multiprocessing, glob
csv.register_dialect("textdialect", delimiter='\t')

def getStats(fn):
	basename = os.path.basename(fn).split('.')[0] + '_stat.txt'
	print "samtools flagstat " + fn + " > " + basename
	os.system("samtools flagstat " + fn + " > " + basename)
	
def main():
	if len(sys.argv)<3: 
		print "Usage: python countMappedReads.py <config file> <output folder>"
		exit()
		
	conf = sys.argv[1]
	ofolder = sys.argv[2]
	if not glob.glob(ofolder + '/'): os.mkdir(ofolder)

	# read list of files that need to be flagstatted
	files = []
	ifile = open(conf, 'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:
		files.append(row[0])
	ifile.close()
	
	# samtools
	os.chdir(ofolder)
	p = multiprocessing.Pool(len(files))
	p.map(getStats, files)
	
	# read into output file
	ofile = open("allstats.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	statfiles = glob.glob("*_stat.txt")
	for fn in statfiles:
		ifile = open(fn, 'r')
		total = ifile.readline().split(' ')[0]
		ifile.readline()
		mapped = ifile.readline().split(' ')[0]
		ifile.close()
		name = fn[:-9]
		writer.writerow([name, total, mapped, float(mapped)/float(total) * 100])
	ofile.close()	
	
if __name__ == "__main__":
	main()
