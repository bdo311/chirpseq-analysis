x=${1%.*}
basename=${x##*/}
fasta=$basename.fa
fastq_trimmed=${basename}_trimmed.fastq
cat $1 | fastq_quality_filter -Q33 -q25 -p80 | fastx_collapser -Q33 > $fasta
perl /seq/scripts/fasta_to_fastq.pl $fasta | fastx_clipper -n -l33 -Q33 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG | fastx_trimmer -f14 -Q33 > $fastq_trimmed
rm -f $fasta
