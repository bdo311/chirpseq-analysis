args = commandArgs(TRUE)
fn = args[1]

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

vals = sapply(vals, function(x)ifelse(is.null(x), 0, x)) #bedgraphs omit positions where read density is 0, so I need to put them back in


## repeats
repeatFn = args[2]
name = args[3]

repeatPos = read.delim(repeatFn, colClasses=c("character", "numeric", "numeric", "numeric"), header=FALSE, row.names=1)
names(repeatPos) = c("len", "start", "stop")
pdf(paste(name, "_repeats.pdf", sep=''), width=15, height=10)
par(mfrow=c(4,5))
for (i in 1:nrow(repeatPos)) {
	start = repeatPos$start[i] + 5
	stop = repeatPos$stop[i] - 5
	maxY = max(vals[start:stop]) * 1.3
	plot(vals[start:stop], type='l', main=rownames(repeatPos)[i], ylim=c(0, maxY), xlab="Position", ylab="Read density")
}
dev.off()


## rRNA
rRNA_start = repeatPos$start[17]
org = args[4]
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
maxY = max(vals[rRNA_start:rRNAend]) * 1.3
plot(vals[rRNA_start:rRNAend], type='l', main="rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")
abline(v=c(start18s, end18s, start5s, end5s, start28s, end28s)-rRNA_start, col="green")
dev.off()


# regions of rRNA
pdf(paste(name, "_rrna_regions.pdf", sep=''), width=8, height=3)
par(mfrow=c(1, 3))

maxY = max(vals[start18s:end18s]) * 1.3
plot(vals[start18s:end18s], type='l', main="18S rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")

maxY = max(vals[start5s:end5s]) * 1.3
plot(vals[start5s:end5s], type='l', main="5S rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")

maxY = max(vals[start28s:end28s]) * 1.3
plot(vals[start28s:end28s], type='l', main="28s rRNA", ylim=c(0, maxY), xlab="Position", ylab="Read density")

dev.off()


