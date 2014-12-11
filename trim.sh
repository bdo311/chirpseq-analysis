x=${1%.*}
basename=${x##*/}
fasta=$basename.fa
fastq_trimmed=${basename}_trimmed.fastq
cat $1 | fastx_trimmer -f7 -Q33 | fastx_clipper -n -l25 -Q33 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG | fastq_quality_filter -Q33 -q25 -p80  > $fastq_trimmed
rm -f $fasta