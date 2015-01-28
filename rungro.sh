#!/bin/sh
# rungro.sh

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
parallel "bowtie2 -p 4 -x $repeat_index -k 1 {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $fastq
parallel "bowtie2 -p 4 -x $genome_index -k 1 {.}_genome.fastq -S {.}_genome.sam" ::: $fastq

#2. samtools - sorting and removing duplicates
parallel "cat {.}_genome.sam | samtools view -Suo - - | samtools sort - {.}_genome_sorted" ::: $fastq
parallel "cat {.}_repeat.sam | samtools view -Suo - - | samtools sort - {.}_repeat_sorted" ::: $fastq
parallel "samtools rmdup -s {.}_genome_sorted.bam {.}_genome_sorted_rmdup.bam" ::: $fastq
parallel "samtools rmdup -s {.}_repeat_sorted.bam {.}_repeat_sorted_rmdup.bam" ::: $fastq

#3. make bedgraph
parallel "bedtools genomecov -ibam {.}_repeat_sorted_rmdup.bam -bg > {.}_repeat.bedGraph; norm_bedGraph.pl {.}_repeat.bedGraph {.}_repeat_norm.bedGraph" ::: $fastq
parallel "bedtools genomecov -ibam {.}_genome_sorted_rmdup.bam -bg > {.}_genome.bedGraph; norm_bedGraph.pl {.}_genome.bedGraph {.}_genome_norm.bedGraph" ::: $fastq
bedGraphToBigWig ${name}_genome_norm.bedGraph $sizes ${name}_genome_norm.bw

#4. get stats for everything
parallel "samtools flagstat {.}_genome_sorted.bam > {.}_genome_sorted_stats.txt" ::: $fastq
parallel "samtools flagstat {.}_repeat_sorted.bam > {.}_repeat_sorted_stats.txt" ::: $fastq

#5. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
parallel "python $scaleScript {.}_repeat_sorted_stats.txt {.}_repeat_norm.bedGraph {.}_repeat_scaled.bedGraph" ::: $fastq
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
Rscript $script $fastq $repeat_pos $name $org

#6. remove sam and bam files
parallel "rm -f {.}_genome_sorted.bam {.}_genome.sam {.}_repeat_sorted.bam {.}_repeat.sam" ::: $fastq
parallel "rm -f {.}_genome.bedGraph" ::: $fastq
parallel "rm -f {.}_repeat_norm.bedGraph {.}_repeat.bedGraph" ::: $fastq

# # exit