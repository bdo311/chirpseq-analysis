#!/bin/sh
# rungro_rnorm.sh

while getopts ":cmh" flag; do
case "$flag" in
    m) bam_only=1
       ;;
	c) compressed=1
	   ;;
	h) echo ""
	   echo "Usage: $0 [-m] [-c] <organism> <name> <reads>"
	   echo ""
	   echo "    -m        Stop after mapping to BAM file"
	   echo "    -c        Input FASTQ is compressed"
	   echo "    -h        Help"
	   echo "    organism  \"mouse\" or \"human\" only"
	   echo "    name      Prefix of output files"
	   echo "    reads     FASTQ or FASTQ.GZ (must specify -c)"  
	   echo ""
	   exit 1
	    ;;
    \?)
       echo "Invalid option: -$OPTARG" >&2
       ;;
  esac
done

shift $((OPTIND-1))

# parameters
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 [-m] [-c] <organism> <name> <reads>"
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
	ets="Mm_CLIP_custom_repeatIndex:2747-6754"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/hg19"
	sizes="/seq/chromosome/hg19/hg19.sizes"
	ets="Repeat_regions:3065-6722"
else
	echo "Organism must be human or mouse."
	exit 1
fi


#1. bowtie
if [ $compressed == 1 ]; then
	prog="zcat"
else 
	prog="cat"
fi

($prog $fastq | bowtie2 -p 8 -x $repeat_index - --un ${name}_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_repeat_sorted) 2> ${name}_bowtie_repeat.err
(cat ${name}_genome.fastq | bowtie2 -p 8 -x $genome_index - | \
samtools view -Suo - - | samtools sort - ${name}_genome_sorted) 2> ${name}_bowtie_genome.err

#2. samtools
samtools index ${name}_genome_sorted.bam
samtools index ${name}_repeat_sorted.bam
samtools flagstat ${name}_genome_sorted.bam > ${name}_genome_sorted_stats.txt
samtools flagstat ${name}_repeat_sorted.bam > ${name}_repeat_sorted_stats.txt
rm -f ${name}_genome.fastq

if [ $bam_only == 1 ]; then
	exit 1
fi

#3: get normalization factor based on number of reads mapping to the 5'ETS
norm_factor="$(samtools view ${name}_repeat_sorted.bam $ets | wc -l)"

#4. bamToBed and split by strand
bedtools bamtobed -i ${name}_genome_sorted.bam | \
awk -F '\t' "BEGIN {OFS=\"\t\"} {if (\$6==\"-\") print \$1,\$2,\$3,\$4,\$5,\$6;}" > \
${name}_genome_sorted_minus.bed
bedtools bamtobed -i ${name}_genome_sorted.bam | \
awk -F '\t' "BEGIN {OFS=\"\t\"} {if (\$6==\"+\") print \$1,\$2,\$3,\$4,\$5,\$6;}" > \
${name}_genome_sorted_plus.bed

#5. output files
bedtools genomecov -i ${name}_genome_sorted_plus.bed -g $sizes -bg > ${name}_genome_plus.bedGraph
gawk -F "\t" -v x=$norm_factor 'BEGIN {OFS="\t"} {print $1,$2,$3,$4 * 1000000 / x}' \
${name}_genome_plus.bedGraph > ${name}_genome_plus_rnorm.bedGraph
norm_bedGraph.pl ${name}_genome_plus.bedGraph ${name}_genome_plus_norm.bedGraph
bedGraphToBigWig ${name}_genome_plus_rnorm.bedGraph $sizes ${name}_genome_plus_rnorm.bw
bedGraphToBigWig ${name}_genome_plus_norm.bedGraph $sizes ${name}_genome_plus_norm.bw

bedtools genomecov -i ${name}_genome_sorted_minus.bed -g $sizes -bg > ${name}_genome_minus.bedGraph
gawk -F "\t" -v x=$norm_factor 'BEGIN {OFS="\t"} {print $1,$2,$3,$4 * 1000000 / x}' \
${name}_genome_minus.bedGraph > ${name}_genome_minus_rnorm.bedGraph
bedGraphToBigWig ${name}_genome_minus_rnorm.bedGraph $sizes ${name}_genome_minus_rnorm.bw
norm_bedGraph.pl ${name}_genome_minus.bedGraph ${name}_genome_minus_norm.bedGraph
bedGraphToBigWig ${name}_genome_minus_norm.bedGraph $sizes ${name}_genome_minus_norm.bw


#6. get stats for everything
samtools flagstat ${name}_genome_sorted.bam > ${name}_genome_sorted_stats.txt
samtools flagstat ${name}_repeat_sorted.bam > ${name}_repeat_sorted_stats.txt

#6. get plots for repeats; this currently doesn't work
bedtools genomecov -ibam ${name}_repeat_sorted.bam -bg > ${name}_repeat.bedGraph
gawk -F "\t" -v x=$norm_factor 'BEGIN {OFS="\t"} {print $1,$2,$3,$4 * 1000000 / x}' \
${name}_repeat.bedGraph > ${name}_repeat_rnorm.bedGraph
# script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
# Rscript $script $fastq $repeat_pos $name $org

#7. remove sam and bam files
rm -f ${name}_genome_plus.bedGraph ${name}_genome_minus.bedGraph
rm -f ${name}_repeat.bedGraph

# # exit