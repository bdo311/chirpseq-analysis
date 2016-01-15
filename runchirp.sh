#!/bin/sh
# runchirp.sh

while getopts ":cmh" flag; do
case "$flag" in
    m) bam_only=1
       ;;
	c) compressed=1
	   ;;
	h) echo ""
	   echo "Usage: $0 [-m] [-c] <organism> <name> <even fastq> <odd fastq> <removal bed file>"
	   echo ""
	   echo "    -m        Stop after mapping to BAM file"
	   echo "    -c        Input FASTQ is compressed"
	   echo "    -h        Help"
	   echo "    organism  \"mouse\" or \"human\" only"
	   echo "    name      Prefix of output files"
	   echo "    fastqs     FASTQ or FASTQ.GZ (must specify -c)"  
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
if [ "$#" -ne 5 ]; then
  echo "$#"
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
	org_macs="mm"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/hg19"
	sizes="/seq/chromosome/hg19/hg19.sizes"
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

($prog $even | bowtie2 -p 8 -x $repeat_index - --un ${name}_even_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_even_repeat_sorted) 2> ${name}_even_bowtie_repeat.err
(cat ${name}_even_genome.fastq | bowtie2 -p 8 -x $genome_index - | \
samtools view -Suo - - | samtools sort - ${name}_even_genome_sorted) 2> ${name}_even_bowtie_genome.err
($prog $odd | bowtie2 -p 8 -x $repeat_index - --un ${name}_odd_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_odd_repeat_sorted) 2> ${name}_odd_bowtie_repeat.err
(cat ${name}_odd_genome.fastq | bowtie2 -p 8 -x $genome_index - | \
samtools view -Suo - - | samtools sort - ${name}_odd_genome_sorted) 2> ${name}_odd_bowtie_genome.err

#2. samtools index/stats
files="${name}_even_genome_sorted.bam ${name}_even_repeat_sorted.bam ${name}_odd_genome_sorted.bam ${name}_odd_repeat_sorted.bam"
parallel "samtools index {}" ::: $files
parallel "samtools flagstat {} > {.}.txt" ::: $files
rm -f ${name}_odd_genome.fastq ${name}_even_genome.fastq

if [ $bam_only == 1 ]; then
	exit 1
fi

#3. getting shifted bedgraphs for reads mapped to genome
python /seq/macs/bin/macs14 -t ${name}_even_genome_sorted.bam -f BAM -n ${name}_even -g $org_macs -B -S
python /seq/macs/bin/macs14 -t ${name}_odd_genome_sorted.bam -f BAM -n ${name}_odd -g $org_macs -B -S
mv ${name}_even_MACS_bedGraph/treat/*.gz ${name}_even_genome_s.bedGraph.gz; gunzip ${name}_even_genome_s.bedGraph.gz
mv ${name}_odd_MACS_bedGraph/treat/*.gz ${name}_odd_genome_s.bedGraph.gz; gunzip ${name}_odd_genome_s.bedGraph.gz
rm -rf ${name}_even_MACS_bedGraph/ ${name}_odd_MACS_bedGraph/

#4. remove RNA component, and normalize
remove_program="/home/raflynn/Scripts/chirpseq_analysis/removeInterval.py"
norm_program="/home/raflynn/Scripts/chirpseq_analysis/normalizeBedgraph.py"
python $remove_program ${name}_odd_genome_s.bedGraph $remove ${name}_odd_genome_sr.bedGraph
python $norm_program ${name}_odd_genome_sr.bedGraph $sizes ${name}_odd_genome_sr_norm.bedGraph
python $remove_program ${name}_even_genome_s.bedGraph $remove ${name}_even_genome_sr.bedGraph
python $norm_program ${name}_even_genome_sr.bedGraph $sizes ${name}_even_genome_sr_norm.bedGraph

# #5. merge genome bedgraph and make bigwig
merge_program="/home/raflynn/Scripts/chirpseq_analysis/takeLower.py"
twofiles="${name}_even_genome_sr_norm.bedGraph ${name}_odd_genome_sr_norm.bedGraph"
bedtools unionbedg -i $twofiles > ${name}_genome_merged_twocol.bedGraph
python $merge_program ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
norm_bedGraph.pl ${name}_genome_merged.bedGraph ${name}_genome_merged_norm.bedGraph
bedGraphToBigWig ${name}_genome_merged_norm.bedGraph $sizes ${name}_genome_merged_norm.bw

#6. make bedgraphs for repeat index; not merging repeat bedgraph
parallel "bedtools genomecov -ibam {.}_repeat_sorted_rmdup.bam -bg > {.}_repeat.bedGraph; norm_bedGraph.pl {.}_repeat.bedGraph {.}_repeat_norm.bedGraph" ::: $even $odd

#7. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
parallel "python $scaleScript {.}_repeat_sorted_rmdup_stats.txt {.}_repeat_norm.bedGraph {.}_repeat_scaled.bedGraph" ::: $even $odd
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
twofiles="${even%%.*}_repeat_scaled.bedGraph ${odd%%.*}_repeat_scaled.bedGraph"
Rscript $script $twofiles $repeat_pos $name $org

#9. remove sam and bam files
rm -f ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
rm -f {.}_genome_shifted.bedGraph {.}_genome_shifted_removed.bedGraph {.}_genome_shifted_removed_norm.bedGraph" ::: $even $odd
parallel "rm -f {.}_repeat_norm.bedGraph {.}_repeat.bedGraph {.}_genome.fastq" ::: $even $odd

# # exit