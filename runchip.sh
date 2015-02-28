#!/bin/sh
# runchip.sh

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

#3. getting shifted bedgraphs for reads mapped to genome
parallel "python /seq/macs/bin/macs14 -t {.}_genome_sorted_rmdup.bam -f BAM -n {.}_genome_shifted -g mm -B -S" ::: $fastq
#parallel "mv {.}_genome_shifted_MACS_bedGraph/treat/*.gz {.}_genome_shifted.bedGraph.gz; rm -rf {.}_genome_shifted_MACS_bedGraph/; gunzip *.gz" ::: $even $odd
parallel "mv {.}_genome_shifted_MACS_bedGraph/treat/*.gz {.}_genome_shifted.bedGraph.gz; gunzip {.}_genome_shifted.bedGraph.gz; rm -rf {.}_genome_shifted_MACS_bedGraph/" ::: $fastq

parallel "norm_bedGraph.pl {.}_genome_shifted.bedGraph ${name}_genome_shifted_norm.bedGraph" ::: $fastq
bedGraphToBigWig ${name}_genome_shifted_norm.bedGraph $sizes ${name}_genome_shifted_norm.bw -clip

#4. make bedgraphs for repeat index
parallel "bedtools genomecov -ibam {.}_repeat_sorted_rmdup.bam -bg > {.}_repeat.bedGraph; norm_bedGraph.pl {.}_repeat.bedGraph {.}_repeat_norm.bedGraph" ::: $fastq

#5. get stats for everything
parallel "samtools flagstat {.}_genome_sorted.bam > {.}_genome_sorted_stats.txt" ::: $fastq
parallel "samtools flagstat {.}_repeat_sorted.bam > {.}_repeat_sorted_stats.txt" ::: $fastq
parallel "samtools flagstat {.}_genome_sorted_rmdup.bam > {.}_genome_sorted_rmdup_stats.txt" ::: $fastq
parallel "samtools flagstat {.}_repeat_sorted_rmdup.bam > {.}_repeat_sorted_rmdup_stats.txt" ::: $fastq

#6. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
parallel "python $scaleScript {.}_repeat_sorted_rmdup_stats.txt {.}_repeat_norm.bedGraph {.}_repeat_scaled.bedGraph" ::: $fastq
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
Rscript $script ${fastq%%.*}_repeat_scaled.bedGraph $repeat_pos $name $org

#7. remove sam and bam files
parallel "rm -f {.}_genome_sorted.bam {.}_genome.sam {.}_repeat_sorted.bam {.}_repeat.sam" ::: $fastq
parallel "rm -f {.}_genome_shifted.bedGraph " ::: $fastq
parallel "rm -f {.}_repeat_norm.bedGraph {.}_repeat.bedGraph {.}_genome.fastq" ::: $fastq
# exit