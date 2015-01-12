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
parallel "bowtie2 -p 4 -x $repeat_index -k 1 -U {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $even $odd
parallel "bowtie2 -p 4 -x $genome_index -k 1 -U {} -S {.}.sam" ::: *_genome.fastq

#2. samtools - sorting and removing duplicates
parallel "cat {} | samtools view -Suo - - | samtools sort - {.}_sorted" ::: *.sam
parallel "samtools rmdup -s {} {.}_rmdup.bam" ::: *_sorted.bam

#3. getting shifted bedgraphs for reads mapped to genome
parallel "python /seq/macs/bin/macs14 -t {.}_sorted_rmdup.bam -f BAM -n {.}_shifted -g mm -B -S" ::: *_genome.sam
parallel "mv {.}_shifted_MACS_bedGraph/treat/*.gz {.}_shifted.bedGraph.gz; rm -rf {.}_shifted_MACS_bedGraph/; gunzip *.gz" ::: *_genome.sam

#4. remove RNA component, and normalize
remove_program="/home/raflynn/Scripts/chirpseq_analysis/removeInterval.py"
norm_program="/home/raflynn/Scripts/chirpseq_analysis/normalizeBedgraph.py"
parallel "python $remove_program {} $remove {.}_removed.bedGraph; python $norm_program {.}_removed.bedGraph $sizes {.}_removed_norm.bedGraph" ::: *_genome_shifted.bedGraph

# #5. merge genome bedgraph and make bigwig
merge_program="/home/raflynn/Scripts/chirpseq_analysis/takeLower.py"
twofiles=$(echo ${name}*_genome_shifted_removed_norm.bedGraph)
bedtools unionbedg -i $twofiles > ${name}_genome_merged_twocol.bedGraph
python $merge_program ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
python $norm_program ${name}_genome_merged.bedGraph $sizes ${name}_genome_merged_norm.bedGraph
bedGraphToBigWig ${name}_genome_merged_norm.bedGraph $sizes ${name}_genome_merged_norm.bw

#6. make bedgraphs for repeat index; not merging repeat bedgraph
parallel "bedtools genomecov -ibam {.}_rmdup.bam -bg > {.}.bedGraph; norm_bedGraph.pl {.}.bedGraph {.}_norm.bedGraph" ::: *_repeat_sorted.bam

#7. get stats for everything
parallel "samtools flagstat {} > {.}_stats.txt" ::: *.bam

#8. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
parallel "python $scaleScript {.}_stats.txt {.}_norm.bedGraph {.}_scaled.bedGraph" ::: *_repeat_sorted.bedGraph
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
twofiles=$(echo *_scaled.bedGraph)
Rscript $script $twofiles $repeat_pos $name $org

# # exit