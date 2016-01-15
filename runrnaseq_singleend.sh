#!/bin/sh
# runrnaseq_singleend.sh

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
	ets="Mm_CLIP_custom_repeatIndex:2747-6754"
	annot="/home/raflynn/Mitosis/RNAseq/annotations/mm9_ucsc-for-tophat.gtf"
	tag="mm9"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/hg19"
	sizes="/seq/chromosome/hg19/hg19.sizes"
	ets="Repeat_regions:3065-6722"
	annot="/home/raflynn/Mitosis/RNAseq/annotations/hg19_ensembl-for-tophat.gtf"
	tag="hg19"
else
	echo "Organism must be human or mouse."
	exit 1
fi

#1. bowtie and tophat
#bowtie2 -p 4 -x $repeat_index -U $fastq --un ${name}_genome.fastq -S ${name}_repeat.sam
bowtie2 -p 4 -x $genome_index -U ${name}_genome.fastq -S ${name}_genome.sam
#tophat -p 4 -o $name -G $annot $tag ${name}_genome.fastq

#1b. cufflinks
#cd $name
#cufflinks -p 4 accepted_hits.bam -g $annot
#cd ..

#2. samtools - sorting
cat ${name}_genome.sam | samtools view -Suo - - | samtools sort - ${name}_genome_sorted
cat ${name}_repeat.sam | samtools view -Suo - - | samtools sort - ${name}_repeat_sorted

#3: get normalization factor based on number of reads mapping to the 5'ETS
samtools index ${name}_repeat_sorted.bam
norm_factor="$(samtools view ${name}_repeat_sorted.bam $ets | wc -l)"

#4. bamToBed and split by strand
bedtools bamtobed -i ${name}_genome_sorted.bam | \
gawk -F '\t' "BEGIN {OFS=\"\t\"} {if (\$6==\"-\") print \$1,\$2,\$3,\$4,\$5,\$6;}" > \
${name}_genome_sorted_minus.bed
bedtools bamtobed -i ${name}_genome_sorted.bam | \
gawk -F '\t' "BEGIN {OFS=\"\t\"} {if (\$6==\"+\") print \$1,\$2,\$3,\$4,\$5,\$6;}" > \
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
norm_bedGraph.pl ${name}_genome_minus.bedGraph ${name}_genome_minus_norm.bedGraph
bedGraphToBigWig ${name}_genome_minus_rnorm.bedGraph $sizes ${name}_genome_minus_rnorm.bw
bedGraphToBigWig ${name}_genome_minus_norm.bedGraph $sizes ${name}_genome_minus_norm.bw

#5. get stats for everything
samtools flagstat ${name}_genome_sorted.bam > ${name}_genome_sorted_stats.txt
samtools flagstat ${name}_repeat_sorted.bam > ${name}_repeat_sorted_stats.txt

#6. remove sam files
rm -f ${name}_repeat.sam ${name}_genome.sam

# # exit