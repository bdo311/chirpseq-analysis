x=${1%.*}
basename=${x##*/}
fasta=$basename.fa
fastq_trimmed=${basename}_trimmed.fastq
cat $1 | fastx_clipper -n -l20 -Q33 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG | fastx_collapser -Q33 > $fasta 
perl /seq/scripts/fasta_to_fastq.pl $fasta > $fastq_trimmed
rm -f $fasta
