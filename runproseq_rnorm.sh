#!/bin/sh
# runproseq_rnorm.sh

# parameters
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <organism> <name> <R1.fastq> <R2.fastq>"
  exit 1
fi

org=$1
name=$2
r1=$3
r2=$4

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
bowtie2 -p 4 -x $repeat_index -1 $r1 -2 $r2 --un-conc ${name}_genome.fastq -S ${name}_repeat.sam
bowtie2 -p 4 -x $genome_index -1 ${name}_genome.1.fastq -2 ${name}_genome.2.fastq -S ${name}_genome.sam

#2. samtools - sorting
cat ${name}_genome.sam | samtools view -Suo - - | samtools sort - ${name}_genome_sorted
cat ${name}_repeat.sam | samtools view -Suo - - | samtools sort - ${name}_repeat_sorted

#3: get normalization factor based on number of reads mapping to the 5'ETS
samtools index ${name}_repeat_sorted.bam
norm_factor="$(samtools view ${name}_repeat_sorted.bam $ets | wc -l)"

#3. bamToBed. Strands are flipped because R2's that map to the '-' strand are at the end of '+' RNAs,
#	and R2's that map to the '+' strand are at the end of '-' RNAs. Only R2's are used because 
#	the 3' end of a nascent transcript is always the first nt of an R2.
bedtools bamtobed -i ${name}_genome_sorted.bam | grep '/2' | \
awk -F '\t' "BEGIN {OFS=\"\t\"} {if (\$6==\"-\") print \$1,\$3-1,\$3,\$4,\$5,\$6;}" > ${name}_genome_sorted_1nt_plus.bed
bedtools bamtobed -i ${name}_genome_sorted.bam | grep '/2' | \
awk -F '\t' "BEGIN {OFS=\"\t\"} {if (\$6==\"+\") print \$1,\$2,\$2+1,\$4,\$5,\$6;}" > ${name}_genome_sorted_1nt_minus.bed

#4. output files
bedtools genomecov -i ${name}_genome_sorted_1nt_plus.bed -g $sizes -bg > ${name}_genome_plus.bedGraph
gawk -F "\t" -v x=$norm_factor 'BEGIN {OFS="\t"} {print $1,$2,$3,$4 * 1000000 / x}' \
${name}_genome_plus.bedGraph > ${name}_genome_plus_rnorm.bedGraph
bedGraphToBigWig ${name}_genome_plus_rnorm.bedGraph $sizes ${name}_genome_plus_rnorm.bw

# bedtools genomecov -i ${name}_genome_sorted_1nt_minus.bed -g $sizes -bg > ${name}_genome_minus.bedGraph
gawk -F "\t" -v x=$norm_factor 'BEGIN {OFS="\t"} {print $1,$2,$3,$4 * 1000000 / x}' \
${name}_genome_minus.bedGraph > ${name}_genome_minus_rnorm.bedGraph
bedGraphToBigWig ${name}_genome_minus_rnorm.bedGraph $sizes ${name}_genome_minus_rnorm.bw

#5. get stats for everything
samtools flagstat ${name}_genome_sorted.bam > ${name}_genome_sorted_stats.txt
samtools flagstat ${name}_repeat_sorted.bam > ${name}_repeat_sorted_stats.txt
