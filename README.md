chirpseq-analysis
=================

This is a set of scripts that together convert ChIRP-seq (chromatin immunoprecipitation followed by RNA purification and sequencing) fastq files into plots and bedgraphs suitable for downstream processing by [metagene-maker](https://github.com/bdo311/metagene-maker/). 

Usage 
------
`chirpseq.sh <organism> <name> <even fastq> <odd fastq> <removal bed file>`

Example: `chirpseq.sh human A375 even_trimmed.fastq odd_trimmed.fastq human_7sk.bed`

### Required arguments:

position | description
------------------|------------------------------------------------
organism | either 'human' or 'mouse'. These are the only two supported at the moment.
name | User-specified prefix for all output files
even fastq, odd fastq | TRIMMED fastq file
removal bed file | genomic positions of the RNA loci that must be masked

Installation instructions and dependencies
--------------

- Clone this repository by running one of the following:
	- `git clone git@github.com:bdo311/chirpseq-analysis.git` if you use ssh authentication
	- `git clone https://github.com/bdo311/chirpseq-analysis.git` otherwise

- Dependencies
  - Bowtie2 (mapping): http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  - MACS 1.4 (peak finding): http://liulab.dfci.harvard.edu/MACS/
  - Samtools (genome manipulation): http://samtools.sourceforge.net/
  - Bedtools (genome arithmetic): http://bedtools.readthedocs.org/en/latest/
  - Python (helper code): https://www.python.org/download/releases/2.7/
  - R (plots): http://www.r-project.org/

Input
-----

2 FASTQ files (even and odd replicates). These FASTQ files must be trimmed (5' and 3'), quality filtered, and collapsed prior to running ChIRPseq analysis.

Output
-----

- normalized bedgraphs for repeat-mapped and genome-mapped reads
- bw file for genome mapped reads
- plots for ChIRP coverage of repeat RNAs
- peak bed file for genome ChIRP


