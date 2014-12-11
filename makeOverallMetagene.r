# makeOverallMetagene.r
# 9/8/14
# uses TSS/TES/gene body files to draw overall metagenes

# Extrapolate bins to exactly the number we need
extrap = function(x, numBins) {
  num = sum(!is.na(x))
  if (num==numBins) return(x)
  approx(1:num, x[1:num], n=numBins)$y
}

# get matrix for files
getMatrix = function(fn, n) {
	no_col = count.fields(fn, sep = "\t") 
	data = t(read.table(fn, sep='\t', fill=TRUE, col.names=1:max(no_col))) #can't use fread

	colnames(data) = data[1,]
	data = data[8:nrow(data),]
	data = apply(data, c(1,2), as.numeric)
	apply(data, 2, extrap, numBins=n) # coerce to exactly some number of bins
}

getMatrix_new = function(fn, start, stop) {
	gro12c.sense.raw = read.table(fn, header=FALSE, nrows=5)
	classes <- sapply(gro12c.sense.raw, class)
	gro12c.sense.raw = read.table(fn, header=FALSE, colClasses=classes)
	genenames = gro12c.sense.raw[,1]
	gro12c.sense.raw = t(gro12c.sense.raw[,start:stop])
	colnames(gro12c.sense.raw) = genenames
	return(t(gro12c.sense.raw))
}

# # divtx and failedterm
# bigenes = read.delim("/home/raflynn/7SK/ChIRPseq/genes/mm9_bidir_genes.txt", sep='\t')
# bigenes.200 = bigenes[bigenes[,6]-bigenes[,3]<200,]
# bidir = c(as.character(bigenes.200[,1]), as.character(bigenes.200[,4]))

# annot = read.delim("/home/raflynn/7SK/ChIRPseq/genes/annotated_genes.txt", sep='\t', header=TRUE)
# div = as.character(annot[annot$udRNA.lfc!=0,]$gene)
# ft = as.character(annot[annot$termination.defect!=0,]$gene)

# tss = getMatrix_new("/home/raflynn/7SK/ChIRPseq/WT2_new/bins/tss/allchr.txt", 8, 207)
# tes = getMatrix_new("/home/raflynn/7SK/ChIRPseq/WT2_new/bins/tes/allchr.txt", 208, 407)
# load("/home/raflynn/7SK/ChIRPseq/WT2_new/bins/geneBody/data.RData")
# geneBody = data.extrap

# tss.all = apply(tss, 2, mean)
# tes.all = apply(tes,2, mean)
# gb.all = apply(geneBody, 1, mean)
# gb.all = approx(1:1000, gb.all[1:1000], n=400)$y
# meta.all = c(tss.all, gb.all, tes.all)

# tss.div = apply(tss[which(rownames(tss) %in% div),], 2, mean)
# tes.div = apply(tes[which(rownames(tes) %in% div),],2, mean)
# gb.div = apply(geneBody[,which(colnames(geneBody) %in% div)], 1, mean)
# gb.div = approx(1:1000, gb.div[1:1000], n=400)$y
# meta.div = c(tss.div, gb.div, tes.div)

# tss.ft = apply(tss[which(rownames(tss) %in% ft),], 2, mean)
# tes.ft = apply(tes[which(rownames(tes) %in% ft),],2, mean)
# gb.ft = apply(geneBody[,which(colnames(geneBody) %in% ft)], 1, mean)
# gb.ft = approx(1:1000, gb.ft[1:1000], n=400)$y
# meta.ft = c(tss.ft, gb.ft, tes.ft)

# pdf("141013_chirp_overall_div_ft.pdf", width=9, height=4.5)
# plot(meta.all, col="red", type='l', lwd=2, xaxt='n', xlab="Position along gene", ylab="7SK read density per nt", ylim=c(0,2.5))
# lines(meta.div, col="green", lwd=2)
# lines(meta.ft, col="gray", lwd=2)
# axis(1, at=c(0, 100, 200, 600, 700, 800), labels=c(-1000, -500, "TSS", "TES", +500, +1000))
# legend("topleft", legend=c("All TSS", "Div. Txn. TSS", "Fail Term. TSS"), lty=1, lwd=2, col=c("red", "green", "gray"))
# dev.off()

getMeta = function(x, start, stop) {
	fns = Sys.glob(x)
	fn=fns[1]
	data=read.delim(fn, row.names=1, header=TRUE)
	return(data[start:stop,])
}

pdf("141126_rnap_overall_metagene.pdf", width=9, height=4.5)

# print("input")
# tss = getMeta("/home/raflynn/mmchipseq/mes_all/Input/tss/avgraw*.txt", 1, 200)
# tes = getMeta("/home/raflynn/mmchipseq/mes_all/Input/tes/avgraw*.txt", 201, 400)
# geneBody = getMeta("/home/raflynn/mmchipseq/mes_all/Input/geneBody/avgraw*.txt", 1, 1000)

# gb.all = approx(1:1000, geneBody, n=400)$y
# input.all = c(tss, gb.all, tes)

plot(1, type='l', xaxt='n', xlab=NA, ylab="7SK read density/nt", xlim=c(0,800), ylim=c(0, 3), lwd=2)
axis(1, at=c(0, 100, 200, 600, 700, 800), labels=c(-1000, -500, "TSS", "TES", +500, +1000))
	
counter=0
#treatments=c("7SK_Input", "WT2_new", "ActD_new", "DRB_new", "Flavo_new", "JQ11_new", "JQ14_new")
#treatments=c("7SK_Input", "WT2_new", "AKT24", "AKT48", "Etoposide", "THZ1", "Triptolide")
treatments=c("Input", "RNAPII_all", "RNAPII_S2", "RNAPII_S5")
#homeDir = "/home/raflynn/7SK/ChIRPseq/"
homeDir = "/home/raflynn/mmchipseq/"
for(treatment in treatments) {
	print(treatment)
	counter=counter + 1
	
	tss = getMeta(paste(homeDir, treatment, "/bins/tss/avgraw*.txt", sep=''), 1,200)
	tes = getMeta(paste(homeDir, treatment, "/bins/tes/avgraw*.txt", sep=''), 201,400)
	geneBody = getMeta(paste(homeDir, treatment, "/bins/geneBody/avgraw*.txt", sep=''),1,1000)

	gb.all = approx(1:1000, geneBody, n=400)$y
	all = c(tss, gb.all, tes)
	
	#all=all/input.all

	lines(all, type='l', col=counter, lwd=2)
}

legend("topright", legend=treatments, lty=1, lwd=2, col=1:length(treatments))
#abline(h=1,lty=2, col="lightblue")
abline(v=c(200,600), lty=2, col="lightblue")
dev.off()












# all-gro

getMatrix_new = function(fn, start, stop) {
	gro12c.sense.raw = read.table(fn, header=FALSE, nrows=5)
	classes <- sapply(gro12c.sense.raw, class)
	gro12c.sense.raw = read.table(fn, header=FALSE, colClasses=classes)
	genenames = gro12c.sense.raw[,1]
	gro12c.sense.raw = t(gro12c.sense.raw[,start:stop])
	colnames(gro12c.sense.raw) = genenames
	return(t(gro12c.sense.raw))
}

pdf("141124_allgenes_gro6hr_comb_overall_metagene.pdf", width=9, height=6)

plot(1, type='l', xaxt='n', xlab=NA, ylab="GROseq read density per nt", xlim=c(0,800), ylim=c(-5, 7), lwd=2)
axis(1, at=c(0, 100, 200, 600, 700, 800), labels=c(-1000, -500, "TSS", "TES", +500, +1000))

treatments=c("6Ccomb", "65comb", "63comb")

counter=0
for(treatment in treatments) {
	print(treatment)
	counter=counter + 1
	print("tss")
	tss = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tss_sense.txt",sep=''), 8, 207)
	print("tes")
	tes = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tes_sense.txt",sep=''), 208, 407)
	print("gb")
	geneBody = read.delim(paste("/home/raflynn/7SK/GROseq/combined/avgraw_GRO_", treatment, "_geneBody_GRO_", treatment, "_geneBody_sense.txt", sep=''), colClasses=c("character", "numeric"), header=TRUE)[,1]
	tss = apply(tss, 2, mean)
	tes = apply(tes, 2, mean)
	geneBody = approx(1:100, geneBody[1:100], n=400)$y

	all = c(tss, geneBody, tes)
	lines(all, type='l', col=counter, lwd=2)
}

counter = 0
for(treatment in treatments) {
	print(treatment)
	counter=counter + 1
	print("tss")
	tss = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tss_antisense.txt",sep=''), 8, 207)
	print("tes")
	tes = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tes_antisense.txt",sep=''), 208, 407)
	print("gb")
	geneBody = read.delim(paste("/home/raflynn/7SK/GROseq/combined/avgraw_GRO_", treatment, "_geneBody_GRO_", treatment, "_geneBody_antisense.txt", sep=''), colClasses=c("character", "numeric"), header=TRUE)[,1]
	tss = apply(tss, 2, mean)
	tes = apply(tes, 2, mean)
	geneBody = approx(1:100, geneBody[1:100], n=400)$y

	all = c(tss, geneBody, tes)
	all = -all
	lines(all, type='l', col=counter, lwd=2)
}
legend("topright", legend=c("6hr Control", "6hr 5' ASO", "6hr 3' ASO"), lty=1, lwd=2, col=c("black", "red", "green"))
abline(h=0)
abline(v=c(200,600),col="lightblue")
dev.off()
