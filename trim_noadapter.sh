x=${1%.*}
basename=${x##*/}
fastq_trimmed=${basename}_trimmed.fastq
cat $1 | fastx_trimmer -f7 -Q33 | fastq_quality_filter -Q33 -q25 -p80 > $fastq_trimmed
