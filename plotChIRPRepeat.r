args = commandArgs(TRUE)
fns = args[1:(length(args)-3)]

fnToVal = list()
counter = 0
for (fn in fns) {
	print(fn)
	data = read.delim(fn, colClasses=c("character", "numeric", "numeric", "numeric"), header=FALSE)
	vals = list()
	for (row in 1:nrow(data)) {
		start = data[row,2]
		end = data[row,3]
		val = data[row,4]
		while(end > start) {
			vals[[start + 1]] = val
			start = start + 1
		}
	}
	counter = counter + 1
	fnToVal[[counter]] = vals
}



## repeats
repeatFn = args[length(args)-2]
repeatPos = read.delim(repeatFn, colClasses=c("character", "numeric", "numeric", "numeric"), header=FALSE, row.names=1)
names(repeatPos) = c("len", "start", "stop")

lastVal = repeatPos[nrow(repeatPos),ncol(repeatPos)]
for (i in 1:length(fnToVal)) {
	vals = fnToVal[[i]]
	vals[[lastVal+6]]=0
	vals = sapply(vals, function(x)ifelse(is.null(x), 0, x)) #bedgraphs omit positions where read density is 0, so I need to put them back in
	fnToVal[[i]] = vals
}

valDF = data.frame(matrix(unlist(fnToVal),nrow=lastVal+6, byrow=F)) #IS NOT HAVING 46069 OK HERE
colnames(valDF)=fns

# SPECIFY COLUMNS WANTED
colsWanted=1:ncol(valDF)
name = args[length(args)-1]
newValDF = valDF

pdf(paste(name, "_repeats.pdf", sep=''), width=15, height=10)
par(mfrow=c(4,5))
for (i in 1:nrow(repeatPos)) {
	start = repeatPos$start[i] - 5
	stop = repeatPos$stop[i] + 5
	maxY = sort(unlist(newValDF[start:stop,colsWanted]), decreasing=TRUE)[10] * 1.3
	plot(1, type='n', main=rownames(repeatPos)[i], xlim=c(0, stop-start), ylim=c(0, maxY), xlab="Position", ylab="Read density")
	for (j in 1:length(colsWanted)) {
		lines(newValDF[start:stop,colsWanted[j]], col=j)	
	}
}
dev.off()


## rRNA
rRNA_start = repeatPos$start[nrow(repeatPos)]
org = args[length(args)]
if (org == "mouse") {
	start18s=rRNA_start +4007
	end18s=rRNA_start +5876
	start5s=rRNA_start +6877
	end5s=rRNA_start +7033
	start28s=rRNA_start +8123
	end28s=rRNA_start +12836
	rRNAend=rRNA_start +13401
} else {
	start18s=rRNA_start +3657
	end18s=rRNA_start +5527
	start5s=rRNA_start +6623
	end5s=rRNA_start +6779
	start28s=rRNA_start +7935
	end28s=rRNA_start +12969
	rRNAend=rRNA_start +13314
}

# whole rRNA
pdf(paste(name, "_rrna.pdf", sep=''), width=8, height=3)

maxY = sort(unlist(newValDF[rRNA_start:rRNAend, colsWanted]), decreasing=TRUE)[10] * 1.3
plot(1, type='n', main="rRNA", xlim=c(0,rRNAend - rRNA_start), ylim=c(0, maxY), xlab="Position", ylab="Read density")
for (j in 1:length(colsWanted)) {
	lines(newValDF[rRNA_start:rRNAend,colsWanted[j]], col=j)
}
	
abline(v=c(start18s, end18s, start5s, end5s, start28s, end28s)-rRNA_start, col="green")
dev.off()

# # regions of rRNA
# pdf(paste(name, "_rrna_regions.pdf", sep=''), width=8, height=3)
# par(mfrow=c(1, 3))

# maxY = max(vals[start18s:end18s]) * 1.3
# plot(vals[start18s:end18s], type='l', main="18S rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")

# maxY = max(vals[start5s:end5s]) * 1.3
# plot(vals[start5s:end5s], type='l', main="5S rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")

# maxY = max(vals[start28s:end28s]) * 1.3
# plot(vals[start28s:end28s], type='l', main="28s rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")

# dev.off()


