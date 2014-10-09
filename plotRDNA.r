fn="input_rdna_sorted.bedGraph"
fn="even_rdna_sorted.bedGraph"
fn="odd_rdna_sorted.bedGraph"
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

ivals = sapply(vals, function(x)ifelse(is.null(x), 0, x)) 
evals = sapply(vals, function(x)ifelse(is.null(x), 0, x))
ovals = sapply(vals, function(x)ifelse(is.null(x), 0, x))

pdf("rdna.pdf", width=6, height=6)
par(mfrow=c(3,1))
plot(ovals[1:13314], type='l', col="red", ylim=c(0,1000), main="Odd 7SK ChIRP", xlab="Position along rRNA", ylab="Read density")
rect(3657,-100,5527,1500,col=rgb(1,0,0,0.3))
rect(6623,-100,6779,1500,col=rgb(0,1,0,0.3))
rect(7935,-100,12969,1500,col=rgb(0,0,1,0.3))

plot(evals[1:13314], type='l', col="blue", ylim=c(0,1000), main="Even 7SK ChIRP", xlab="Position along rRNA", ylab="Read density")
rect(3657,-100,5527,1500,col=rgb(1,0,0,0.3))
rect(6623,-100,6779,1500,col=rgb(0,1,0,0.3))
rect(7935,-100,12969,1500,col=rgb(0,0,1,0.3))

plot(ivals[1:13314], type='l', col="black", ylim=c(0,1000), main="Input 7SK ChIRP", xlab="Position along rRNA", ylab="Read density")
rect(3657,-100,5527,1500,col=rgb(1,0,0,0.3))
rect(6623,-100,6779,1500,col=rgb(0,1,0,0.3))
rect(7935,-100,12969,1500,col=rgb(0,0,1,0.3))

dev.off()