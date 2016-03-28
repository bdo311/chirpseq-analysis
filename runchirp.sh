#!/bin/sh
# runchirp.sh

bam_only=0
compressed=0
to_remove=0

while getopts ":cmhr:" flag; do
case "$flag" in
    m) bam_only=1
       ;;
	c) compressed=1
	   ;;
	r) removal_file="$OPTARG"
	   to_remove=1
	   ;;
	h) echo ""
	   echo "Usage: sh $0 [-m] [-c] [-r removal bed file] <organism> <name> <STAR index> <even fastq> <odd fastq> <input fastq>"
	   echo ""
	   echo "    -m                  Stop after mapping to BAM file"
	   echo "    -c                  Input FASTQ is compressed"
	   echo "    -h                  Help"
	   echo "    -r REMOVAL_FILE     File to remove ChIRPseq reads from"
	   echo "    organism            \"mouse\" or \"human\" only"
	   echo "    name                Prefix of output files"
	   echo "    STAR index          Path to STAR index"
	   echo "    fastqs              FASTQ or FASTQ.GZ (must specify -c if .gz)"  
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
if [ "$#" -ne 6 ]; then
  echo "Usage: sh $0 [-m] [-c] [-r removal bed file] <organism> <name> <STAR index> <even fastq> <odd fastq> <input fastq>"
  exit 1
fi

org=$1
name=$2
genome_index=$3
even=$4
odd=$5
input=$6

if [ $org == "mouse" ]; then
	repeat_pos="~/Scripts/repeat_index/mm9/Mm_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/mm9/rep_spaced"
	sizes="/seq/chromosome/mm9/mm9.sizes"
	org_macs="mm"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	sizes="/seq/chromosome/hg19/hg19.sizes"
	org_macs="hs"
else
	echo "Organism must be human or mouse."
	exit 1
fi

#1. bowtie
if [ $compressed -eq 1 ]; then
	prog="zcat"
else 
	prog="cat"
fi

# bowtie to map against repeats
($prog $even | bowtie2 -p 8 -x $repeat_index - --un ${name}_even_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_even_repeat_sorted) 2> ${name}_even_bowtie_repeat.err
($prog $odd | bowtie2 -p 8 -x $repeat_index - --un ${name}_odd_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_odd_repeat_sorted) 2> ${name}_odd_bowtie_repeat.err
($prog $input | bowtie2 -p 8 -x $repeat_index - --un ${name}_input_genome.fastq | \
samtools view -Suo - - | samtools sort - ${name}_input_repeat_sorted) 2> ${name}_input_bowtie_repeat.err

# STAR to map against genome
STAR --genomeDir $genome_index --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn ${name}_even_genome.fastq \
--outFileNamePrefix ${name}_even_genome_ --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.08
STAR --genomeDir $genome_index --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn ${name}_odd_genome.fastq \
--outFileNamePrefix ${name}_odd_genome_ --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.08
STAR --genomeDir $genome_index --runThreadN 8 --genomeLoad LoadAndKeep --readFilesIn ${name}_input_genome.fastq \
--outFileNamePrefix ${name}_input_genome_ --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.08

# using samtools because STAR can't output to command line
cat ${name}_even_genome_Aligned.out.sam | samtools view -Su -q 255 - -o - | samtools sort - ${name}_even_genome_sorted
cat ${name}_odd_genome_Aligned.out.sam | samtools view -Su -q 255 - -o - | samtools sort - ${name}_odd_genome_sorted
cat ${name}_input_genome_Aligned.out.sam | samtools view -Su -q 255 - -o - | samtools sort - ${name}_input_genome_sorted

# (cat ${name}_even_genome.fastq | bowtie2 -p 8 -x $genome_index - | \
# samtools view -Suo - - | samtools sort - ${name}_even_genome_sorted) 2> ${name}_even_bowtie_genome.err
# (cat ${name}_odd_genome.fastq | bowtie2 -p 8 -x $genome_index - | \
# samtools view -Suo - - | samtools sort - ${name}_odd_genome_sorted) 2> ${name}_odd_bowtie_genome.err

#2. samtools index/stats
files="${name}_even_genome_sorted.bam ${name}_even_repeat_sorted.bam ${name}_odd_genome_sorted.bam \
${name}_odd_repeat_sorted.bam ${name}_input_repeat_sorted.bam ${name}_input_genome_sorted.bam"
parallel "samtools index {}" ::: $files
parallel "samtools flagstat {} > {.}_stats.txt" ::: $files
if [ ! -d "todelete" ]; then
  mkdir todelete
fi
mv ${name}_odd_genome.fastq ${name}_even_genome.fastq ${name}_input_genome.fastq todelete

if [ $bam_only == 1 ]; then
	exit 1
fi

#3. getting shifted bedgraphs for reads mapped to genome
python /seq/macs/bin/macs14 -t ${name}_even_genome_sorted.bam -f BAM -n ${name}_even -g $org_macs -B -S
python /seq/macs/bin/macs14 -t ${name}_odd_genome_sorted.bam -f BAM -n ${name}_odd -g $org_macs -B -S
python /seq/macs/bin/macs14 -t ${name}_input_genome_sorted.bam -f BAM -n ${name}_input -g $org_macs -B -S
mv ${name}_even_MACS_bedGraph/treat/*.gz ${name}_even_genome_s.bedGraph.gz; gunzip ${name}_even_genome_s.bedGraph.gz
mv ${name}_odd_MACS_bedGraph/treat/*.gz ${name}_odd_genome_s.bedGraph.gz; gunzip ${name}_odd_genome_s.bedGraph.gz
mv ${name}_input_MACS_bedGraph/treat/*.gz ${name}_input_genome_s.bedGraph.gz; gunzip ${name}_input_genome_s.bedGraph.gz
rm -rf ${name}_even_MACS_bedGraph/ ${name}_odd_MACS_bedGraph/ ${name}_input_MACS_bedGraph/

#4. remove RNA component, and normalize
remove_program="/home/raflynn/Scripts/chirpseq_analysis/removeInterval.py"
norm_program="/home/raflynn/Scripts/chirpseq_analysis/normalizeBedgraph.py"

if [ $to_remove == 1 ]; then
	python $remove_program ${name}_odd_genome_s.bedGraph $removal_file ${name}_odd_genome_sr.bedGraph
	mv ${name}_odd_genome_sr.bedGraph ${name}_odd_genome_s.bedGraph
	python $norm_program ${name}_odd_genome_s.bedGraph $sizes ${name}_odd_genome_sr_norm.bedGraph
	python $remove_program ${name}_even_genome_s.bedGraph $removal_file ${name}_even_genome_sr.bedGraph
	mv ${name}_even_genome_sr.bedGraph ${name}_even_genome_s.bedGraph
	python $norm_program ${name}_even_genome_s.bedGraph $sizes ${name}_even_genome_sr_norm.bedGraph
	python $norm_program ${name}_input_genome_s.bedGraph $sizes ${name}_input_genome_sr_norm.bedGraph  # no need to remove
else 
	python $norm_program ${name}_odd_genome_s.bedGraph $sizes ${name}_odd_genome_sr_norm.bedGraph
	python $norm_program ${name}_even_genome_s.bedGraph $sizes ${name}_even_genome_sr_norm.bedGraph
	python $norm_program ${name}_input_genome_s.bedGraph $sizes ${name}_input_genome_sr_norm.bedGraph
fi

#5. merge genome bedgraph and make bigwig
merge_program="/home/raflynn/Scripts/chirpseq_analysis/takeLower.py"
twofiles="${name}_even_genome_sr_norm.bedGraph ${name}_odd_genome_sr_norm.bedGraph"
bedtools unionbedg -i $twofiles > ${name}_genome_merged_twocol.bedGraph
python $merge_program ${name}_genome_merged_twocol.bedGraph ${name}_genome_merged.bedGraph
norm_bedGraph.pl ${name}_genome_merged.bedGraph ${name}_genome_merged_norm.bedGraph
bedGraphToBigWig ${name}_genome_merged_norm.bedGraph $sizes ${name}_genome_merged_norm.bw

#6. make bedgraphs for repeat index; not merging repeat bedgraph
bedtools genomecov -ibam ${name}_even_repeat_sorted.bam -bg > ${name}_even_repeat.bedGraph
norm_bedGraph.pl ${name}_even_repeat.bedGraph ${name}_even_repeat_norm.bedGraph
bedtools genomecov -ibam ${name}_odd_repeat_sorted.bam -bg > ${name}_odd_repeat.bedGraph
norm_bedGraph.pl ${name}_odd_repeat.bedGraph ${name}_odd_repeat_norm.bedGraph

#7. get plots for repeats
scaleScript="/home/raflynn/Scripts/chirpseq_analysis/rescaleRepeatBedgraph.py"
parallel "python $scaleScript {}_repeat_sorted_stats.txt {}_repeat_norm.bedGraph {}_repeat_scaled.bedGraph" ::: ${name}_even ${name}_odd
script="/home/raflynn/Scripts/chirpseq_analysis/plotChIRPRepeat.r"
twofiles="${name}_even_repeat_scaled.bedGraph ${name}_odd_repeat_scaled.bedGraph"
Rscript $script $twofiles $repeat_pos $name $org

#8. get a stringent peak list
bedGraph2sam_script="/home/raflynn/Scripts/chirpseq_analysis/bedGraph2sam.pl"
peak_correlation_script="/home/raflynn/Scripts/chirpseq_analysis/peak_correlation.pl"
twofiles="${name}_even_genome_s.bedGraph ${name}_odd_genome_s.bedGraph"  # not normalized
bedtools unionbedg -i $twofiles > ${name}_genome_merged_twocol_unnorm.bedGraph
python $merge_program ${name}_genome_merged_twocol_unnorm.bedGraph ${name}_genome_merged_unnorm.bedGraph
perl $bedGraph2sam_script $sizes ${name}_genome_merged_unnorm.bedGraph ${name}_genome_merged_unnorm.sam
python /seq/macs/bin/macs14 -t ${name}_genome_merged_unnorm.sam -c ${name}_input_genome_Aligned.out.sam -f SAM -n ${name}_genome_merge -g $org_macs --bw=200 -m 10,50
perl $peak_correlation_script ${name}_genome_merge_peaks.xls ${name}_even_genome_sr_norm.bedGraph ${name}_odd_genome_sr_norm.bedGraph ${name}_genome_merged_norm.bedGraph
mv ${name}_genome_merge_peaks.xls.corr.xls ${name}_genome_merged_confidentpeaks.xls

#9. remove sam and bam files
mv *_s.bedGraph *_sr.bedGraph todelete
mv *.sam todelete
mv *_unnorm.bedGraph todelete
mv *_repeat.bedGraph todelete
mv *genome_sr_norm.bedGraph todelete
mv *_Log.out *.progress.out *.tab todelete
mv *.tab todelete

exit 0
