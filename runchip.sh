#!/bin/sh
# runchip_compressed.sh

bam_only=0
compressed=0

while getopts ":cmh" flag; do
case "$flag" in
    m) bam_only=1
       ;;
	c) compressed=1
	   ;;
	h) echo ""
	   echo "Usage: $0 [-m] [-c] <organism> <name> <STAR index> <reads>"
	   echo ""
	   echo "    -m                  Stop after mapping to BAM file"
	   echo "    -c                  Input FASTQ is compressed"
	   echo "    -h                  Help"
	   echo "    organism            \"mouse\" or \"human\" only"
	   echo "    name                Prefix of output files"
	   echo "    STAR index          Path to STAR index"
	   echo "    reads               FASTQ or FASTQ.GZ (must specify -c)"  
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
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 [-m] [-c] <organism> <name> <STAR index> <reads>"
  exit 1
fi

org=$1
name=$2
genome_index=$3
fastq=$4

if [ $org == "mouse" ]; then
	repeat_pos="~/Scripts/repeat_index/mm9/Mm_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/mm9/rep_spaced"
	sizes="/seq/chromosome/mm9/mm9.sizes"
	ets="Mm_CLIP_custom_repeatIndex:2747-6754"
	org_macs="mm"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	sizes="/seq/chromosome/hg19/hg19.sizes"
	ets="Repeat_regions:3065-6722"
	org_macs="hs"
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

# 1. map to repeat and genome

# bowtie to map against repeat index
($prog $fastq | bowtie2 -p 8 -x $repeat_index - --un ${name}_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_repeat_sorted) 2> ${name}_bowtie_repeat.err
# (cat ${name}_genome.fastq | bowtie2 -p 8 -x $genome_index - | \
# samtools view -Suo - - | samtools sort - ${name}_genome_sorted) 2> ${name}_bowtie_genome.err

# STAR to map against genome
STAR --genomeDir $genome_index --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn ${name}_genome.fastq \
--outFileNamePrefix ${name}_genome_ --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.08
cat ${name}_genome_Aligned.out.sam | samtools view -Su -q 255 - -o - | samtools sort - ${name}_genome_sorted

#2. samtools index/stats
samtools index ${name}_genome_sorted.bam
samtools index ${name}_repeat_sorted.bam
samtools flagstat ${name}_genome_sorted.bam > ${name}_genome_sorted_stats.txt
samtools flagstat ${name}_repeat_sorted.bam > ${name}_repeat_sorted_stats.txt
rm -f ${name}_genome.fastq

if [ $bam_only == 1 ]; then
	exit 1
fi

#3. getting shifted bedgraphs for reads mapped to genome
python /seq/macs/bin/macs14 -t ${name}_genome_sorted.bam -f BAM -n ${name}_genome_shifted -g $org_macs -B -S
mv ${name}_genome_shifted_MACS_bedGraph/treat/*.gz ${name}_genome_shifted.bedGraph.gz
gunzip ${name}_genome_shifted.bedGraph.gz
rm -rf ${name}_genome_shifted_MACS_bedGraph/

samtools index ${name}_repeat_sorted.bam
norm_factor="$(samtools view ${name}_repeat_sorted.bam $ets | wc -l)"

gawk -F "\t" -v x=$norm_factor 'BEGIN {OFS="\t"} {print $1,$2,$3,$4 * 1000000 / x}' \
${name}_genome_shifted.bedGraph > ${name}_genome_shifted_rnorm.bedGraph
bedGraphToBigWig ${name}_genome_shifted_rnorm.bedGraph $sizes ${name}_genome_shifted_rnorm.bw -clip
norm_bedGraph.pl ${name}_genome_shifted.bedGraph ${name}_genome_shifted_norm.bedGraph
bedGraphToBigWig ${name}_genome_shifted_norm.bedGraph $sizes ${name}_genome_shifted_norm.bw -clip

#4. make bedgraphs for repeat index
bedtools genomecov -ibam ${name}_repeat_sorted.bam -bg > ${name}_repeat.bedGraph; norm_bedGraph.pl ${name}_repeat.bedGraph ${name}_repeat_norm.bedGraph

#5. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
python $scaleScript ${name}_repeat_sorted_stats.txt ${name}_repeat_norm.bedGraph ${name}_repeat_scaled.bedGraph
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
Rscript $script ${name}_repeat_scaled.bedGraph $repeat_pos $name $org

#6. remove unneeded files
rm -f *_Log.progress.out *Log.out *_SJ.out.tab
rm -f *_genome_shifted.bedGraph
rm -f *.sam

exit 0
