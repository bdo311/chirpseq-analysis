#!/bin/sh
# runchipnexus.sh

# parameters
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <organism> <name> <fastq>"
  exit 1
fi

org=$1
name=$2
fastq=$3

if [ $org == "mouse" ]; then
	repeat_pos="~/Scripts/repeat_index/mm9/Mm_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/mm9/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/mm9"
	sizes="/seq/chromosome/mm9/mm9.sizes"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/hg19"
	sizes="/seq/chromosome/hg19/hg19.sizes"
else
	echo "Organism must be human or mouse."
	exit 1
fi

#1. bowtie
parallel "bowtie2 -p 4 -x $repeat_index {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $fastq
parallel "bowtie2 -p 4 -x $genome_index {.}_genome.fastq -S {.}_genome.sam" ::: $fastq

#2. samtools - sorting and removing duplicates
parallel "cat {.}_genome.sam | samtools view -Suo - - | samtools sort - {.}_genome_sorted" ::: $fastq
parallel "cat {.}_repeat.sam | samtools view -Suo - - | samtools sort - {.}_repeat_sorted" ::: $fastq

#3. make bed file
parallel "bedtools bamtobed -i {.}_genome_sorted.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {if (\$6==\"+\") print \$1,\$2,\$2+1,\$4,\$5,\$6;}' > {.}_genome_sorted_1nt_plus.bed" ::: $fastq
parallel "bedtools bamtobed -i {.}_genome_sorted.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {if (\$6==\"-\") print \$1,\$3-1,\$3,\$4,\$5,\$6;}' > {.}_genome_sorted_1nt_minus.bed" ::: $fastq
parallel "bedtools bamtobed -i {.}_repeat_sorted.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {print \$1,\$2,\$2+1,\$4,\$5,\$6}' > {.}_repeat_sorted_1nt.bed" ::: $fastq

#4. make bedgraphs 
awk 'BEGIN {print "Repeat_regions\t50000"}' > repeat_regions.txt
parallel "bedtools genomecov -i {.}_repeat_sorted_1nt.bed -g repeat_regions.txt -bg > {.}_repeat.bedGraph; norm_bedGraph.pl {.}_repeat.bedGraph {.}_repeat_norm.bedGraph" ::: $fastq
rm -f repeat_regions.txt
parallel "bedtools genomecov -i {.}_genome_sorted_1nt_plus.bed -g $sizes -bg > {.}_genome_plus.bedGraph; norm_bedGraph.pl {.}_genome_plus.bedGraph {.}_genome_plus_norm.bedGraph" ::: $fastq
parallel "bedGraphToBigWig {.}_genome_plus_norm.bedGraph $sizes {.}_genome_plus_norm.bw" ::: $fastq
parallel "bedtools genomecov -i {.}_genome_sorted_1nt_minus.bed -g $sizes -bg > {.}_genome_minus.bedGraph; norm_bedGraph.pl {.}_genome_minus.bedGraph {.}_genome_minus_norm.bedGraph" ::: $fastq
parallel "bedGraphToBigWig {.}_genome_minus_norm.bedGraph $sizes {.}_genome_minus_norm.bw" ::: $fastq

#5. get stats for everything
parallel "samtools flagstat {.}_genome_sorted.bam > {.}_genome_sorted_stats.txt" ::: $fastq
parallel "samtools flagstat {.}_repeat_sorted.bam > {.}_repeat_sorted_stats.txt" ::: $fastq

#7. remove sam and bam files
#parallel "rm -f {.}_genome_sorted.bam {.}_genome.sam {.}_repeat_sorted.bam {.}_repeat.sam" ::: $fastq
#parallel "rm -f {.}_genome_shifted.bedGraph " ::: $fastq
#parallel "rm -f {.}_repeat_norm.bedGraph {.}_repeat.bedGraph {.}_genome.fastq" ::: $fastq
# exit