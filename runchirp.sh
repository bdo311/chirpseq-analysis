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
	repeat_pos="~/Scripts/repeat_index/mm9/Mm_repeatIndex_positions.txt"
	repeat_index="~/Scripts/repeat_index/mm9/rep"
	genome_index="/seq/bowtie2-2.1.0/indexes/mm9"
	sizes="/seq/chromosome/mm9/mm9.sizes"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep"
	genome_index="/seq/bowtie2-2.1.0/indexes/hg19"
	sizes="/seq/chromosome/hg19/hg19.sizes"
else
	echo "Organism must be human or mouse."
	exit 1
fi

#1. bowtie
parallel "bowtie2 -p 4 -x $repeat_index -k 1 -U {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $even $odd
parallel "bowtie2 -p 4 -x $genome_index -k 1 -U {} -S {.}.sam" ::: *_genome.fastq

#2. getting shifted bedgraphs for reads mapped to genome
parallel "macs14 -t {} -n {.}_shifted -g mm -B -S" ::: *_genome.sam
parallel "mv {.}_shifted_MACS_bedGraph/treat/*.gz {.}_shifted.bedGraph.gz; rm -rf {.}_shifted_MACS_bedGraph/; gunzip *.gz" ::: *_genome.sam

# 3. remove RNA component, and normalize
remove_program="/home/raflynn/Scripts/chirpseq_analysis/removeInterval.py"
norm_program="/home/raflynn/Scripts/chirpseq_analysis/normalizeBedgraph.py"
parallel "python $remove_program {} $remove {.}_removed.bedGraph; python $norm_program {.}_removed.bedGraph $sizes {.}_removed_norm.bedGraph" ::: *_genome_shifted.bedGraph

# 4. merge genome bedgraph
merge_program="/home/raflynn/Scripts/chirpseq_analysis/takeLower.py"
twofiles=$(echo *_genome_shifted_removed_norm.bedGraph)
bedtools unionbedg -i $twofiles > genome_merged_twocol.bedGraph
python $merge_program genome_merged_twocol.bedGraph genome_merged.bedGraph
python $norm_program genome_merged.bedGraph $sizes ${name}_genome_merged_norm.bedGraph

# 5. bigwig for genome
bedGraphToBigWig ${name}_genome_merged_norm.bedGraph $sizes ${name}_genome_merged_norm.bw

#6. samtools to make bedgraphs for repeat index
parallel "samtools view -Suo - {} | samtools sort - {.}_sorted" ::: *_repeat.sam
parallel "bedtools genomecov -ibam {} -bg > {.}.bedGraph; norm_bedGraph.pl {.}.bedGraph {.}_norm.bedGraph" ::: *_repeat_sorted.bam

# 7. merge repeat bedgraph
twofiles=$(echo *_repeat_sorted_norm.bedGraph)
bedtools unionbedg -i $twofiles > repeat_merged_twocol.bedGraph
python $merge_program repeat_merged_twocol.bedGraph repeat_merged.bedGraph
norm_bedGraph.pl repeat_merged.bedGraph  ${name}_repeat_merged_norm.bedGraph

#8. get plot for repeats
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
Rscript $script ${name}_repeat_merged_norm.bedGraph $repeat_pos $name $org

# # exit