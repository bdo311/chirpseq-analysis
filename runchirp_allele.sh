#!/bin/sh
# runchirp_allele.sh

# parameters
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <organism> <name> <even bam> <odd bam> <removal bed file>"
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

#3. getting shifted bedgraphs for reads mapped to genome
if [ $org == "mouse" ]; then
	macs2 callpeak -f BAMPE -t $even -n ${name}_even -B -g mm
	macs2 callpeak -f BAMPE -t $odd -n ${name}_odd -B -g mm
elif [ $org == "human" ]; then
	macs2 callpeak -f BAMPE -t $even -n ${name}_even -B
	macs2 callpeak -f BAMPE -t $odd -n ${name}_odd -B
else
	echo "Organism must be human or mouse."
	exit 1	
fi
mv ${name}_even_treat_pileup.bdg ${name}_even_genome_shifted.bedGraph
mv ${name}_odd_treat_pileup.bdg ${name}_odd_genome_shifted.bedGraph
rm -f ${name}_even_control_lambda.bdg ${name}_odd_control_lambda.bdg 

#4. remove RNA component, and normalize
remove_program="/home/raflynn/Scripts/chirpseq_analysis/removeInterval.py"
norm_program="/home/raflynn/Scripts/chirpseq_analysis/normalizeBedgraph.py"
python $remove_program ${name}_even_genome_shifted.bedGraph $remove ${name}_even_genome_shifted_removed.bedGraph
python $norm_program ${name}_even_genome_shifted_removed.bedGraph $sizes ${name}_even_genome_shifted_removed_norm.bedGraph
python $remove_program ${name}_odd_genome_shifted.bedGraph $remove ${name}_odd_genome_shifted_removed.bedGraph
python $norm_program ${name}_odd_genome_shifted_removed.bedGraph $sizes ${name}_odd_genome_shifted_removed_norm.bedGraph
sort -k1,1 -k2,2n ${name}_even_genome_shifted_removed_norm.bedGraph > ${name}_even_genome_shifted_removed_norm_st.bedGraph
sort -k1,1 -k2,2n ${name}_odd_genome_shifted_removed_norm.bedGraph > ${name}_odd_genome_shifted_removed_norm_st.bedGraph

# #5. merge genome bedgraph and make bigwig
merge_program="/home/raflynn/Scripts/chirpseq_analysis/takeLower.py"
twofiles="${name}_even_genome_shifted_removed_norm_st.bedGraph ${name}_odd_genome_shifted_removed_norm_st.bedGraph"
bedtools unionbedg -i $twofiles > ${name}_genome_merged_twocol.bedGraph
python $merge_program ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
norm_bedGraph.pl ${name}_genome_merged.bedGraph ${name}_genome_merged_norm.bedGraph
bedGraphToBigWig ${name}_genome_merged_norm.bedGraph $sizes ${name}_genome_merged_norm.bw

# #6. make bedgraphs for repeat index; not merging repeat bedgraph
# parallel "bedtools genomecov -ibam {.}_repeat_sorted_rmdup.bam -bg > {.}_repeat.bedGraph; norm_bedGraph.pl {.}_repeat.bedGraph {.}_repeat_norm.bedGraph" ::: $even $odd

# #7. get stats for everything
# parallel "samtools flagstat {} > {.}_genome_sorted_stats.txt" ::: $even $odd

# #8. get plots for repeats
# scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
# parallel "python $scaleScript {.}_repeat_sorted_rmdup_stats.txt {.}_repeat_norm.bedGraph {.}_repeat_scaled.bedGraph" ::: $even $odd
# script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
# twofiles="${even%%.*}_repeat_scaled.bedGraph ${odd%%.*}_repeat_scaled.bedGraph"
# Rscript $script $twofiles $repeat_pos $name $org

#9. remove sam and bam files
rm -f ${name}_genome_merged_twocol.bedGraph 
rm -f ${name}_even_genome_shifted.bedGraph ${name}_even_genome_shifted_removed.bedGraph ${name}_even_genome_shifted_removed_norm.bedGraph ${name}_odd_genome_shifted_removed_norm_st.bedGraph
rm -f ${name}_odd_genome_shifted.bedGraph ${name}_odd_genome_shifted_removed.bedGraph ${name}_odd_genome_shifted_removed_norm.bedGraph ${name}_odd_genome_shifted_removed_norm_st.bedGraph

# # exit