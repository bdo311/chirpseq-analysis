#!/bin/sh
# runchirp.sh

# parameters
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <organism> <name> <even fastq> <odd fastq> <removal bed file>"
  exit 1
fi

org=$1
name=$2
even=$3
odd=$4
remove=$5

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
parallel "bowtie2 -p 4 -x $repeat_index {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $even $odd
parallel "bowtie2 -p 4 -x $genome_index {.}_genome.fastq -S {.}_genome.sam" ::: $even $odd

#2. samtools - sorting and removing duplicates
parallel "cat {.}_genome.sam | samtools view -Suo - - | samtools sort - {.}_genome_sorted" ::: $even $odd
parallel "cat {.}_repeat.sam | samtools view -Suo - - | samtools sort - {.}_repeat_sorted" ::: $even $odd
parallel "samtools rmdup -s {.}_genome_sorted.bam {.}_genome_sorted_rmdup.bam" ::: $even $odd
parallel "samtools rmdup -s {.}_repeat_sorted.bam {.}_repeat_sorted_rmdup.bam" ::: $even $odd

#3. getting shifted bedgraphs for reads mapped to genome
parallel "python /seq/macs/bin/macs14 -t {.}_genome_sorted_rmdup.bam -f BAM -n {.}_genome_shifted -g mm -B -S" ::: $even $odd
#parallel "mv {.}_genome_shifted_MACS_bedGraph/treat/*.gz {.}_genome_shifted.bedGraph.gz; rm -rf {.}_genome_shifted_MACS_bedGraph/; gunzip *.gz" ::: $even $odd
parallel "mv {.}_genome_shifted_MACS_bedGraph/treat/*.gz {.}_genome_shifted.bedGraph.gz; gunzip {.}_genome_shifted.bedGraph.gz" ::: $even $odd

#4. remove RNA component, and normalize
remove_program="/home/raflynn/Scripts/chirpseq_analysis/removeInterval.py"
norm_program="/home/raflynn/Scripts/chirpseq_analysis/normalizeBedgraph.py"
parallel "python $remove_program {.}_genome_shifted.bedGraph $remove {.}_genome_shifted_removed.bedGraph; python $norm_program {.}_genome_shifted_removed.bedGraph $sizes {.}_genome_shifted_removed_norm.bedGraph" ::: $even $odd

# #5. merge genome bedgraph and make bigwig
merge_program="/home/raflynn/Scripts/chirpseq_analysis/takeLower.py"
twofiles="${even%%.*}_genome_shifted_removed_norm.bedGraph ${odd%%.*}_genome_shifted_removed_norm.bedGraph"
bedtools unionbedg -i $twofiles > ${name}_genome_merged_twocol.bedGraph
python $merge_program ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
norm_bedGraph.pl ${name}_genome_merged.bedGraph ${name}_genome_merged_norm.bedGraph
bedGraphToBigWig ${name}_genome_merged_norm.bedGraph $sizes ${name}_genome_merged_norm.bw

#6. make bedgraphs for repeat index; not merging repeat bedgraph
parallel "bedtools genomecov -ibam {.}_repeat_sorted_rmdup.bam -bg > {.}_repeat.bedGraph; norm_bedGraph.pl {.}_repeat.bedGraph {.}_repeat_norm.bedGraph" ::: $even $odd

#7. get stats for everything
parallel "samtools flagstat {.}_genome_sorted.bam > {.}_genome_sorted_stats.txt" ::: $even $odd
parallel "samtools flagstat {.}_repeat_sorted.bam > {.}_repeat_sorted_stats.txt" ::: $even $odd
parallel "samtools flagstat {.}_genome_sorted_rmdup.bam > {.}_genome_sorted_rmdup_stats.txt" ::: $even $odd
parallel "samtools flagstat {.}_repeat_sorted_rmdup.bam > {.}_repeat_sorted_rmdup_stats.txt" ::: $even $odd

#8. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
parallel "python $scaleScript {.}_repeat_sorted_rmdup_stats.txt {.}_repeat_norm.bedGraph {.}_repeat_scaled.bedGraph" ::: $even $odd
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
twofiles="${even%%.*}_repeat_scaled.bedGraph ${odd%%.*}_repeat_scaled.bedGraph"
Rscript $script $twofiles $repeat_pos $name $org

#9. remove sam and bam files
parallel "rm -f {.}_genome_sorted.bam {.}_genome.sam {.}_repeat_sorted.bam {.}_repeat.sam" ::: $even $odd
rm -f ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
parallel "rm -f {.}_genome_shifted.bedGraph {.}_genome_shifted_removed.bedGraph {.}_genome_shifted_removed_norm.bedGraph" ::: $even $odd
parallel "rm -f {.}_repeat_norm.bedGraph {.}_repeat.bedGraph {.}_genome.fastq" ::: $even $odd

# # exit