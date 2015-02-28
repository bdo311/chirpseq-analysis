# averageBedgraphs.sh
# 2/27/15
# usage: sh averageBedgraphs.sh bg1 bg2 bg3 ... output.bedgraph

array=( $@ )
len=${#array[@]}
_output=${array[$len-1]}
_args=${array[@]:0:$len-1}
n=$RANDOM

merge_program="/home/raflynn/Scripts/chirpseq_analysis/makeAverageBedgraph.py"
bedtools unionbedg -i $_args > ${n}_genome_merged_col.bedGraph
python $merge_program ${n}_genome_merged_col.bedGraph ${n}_genome_merged.bedGraph
rm -f ${n}_genome_merged_col.bedGraph
norm_bedGraph.pl ${n}_genome_merged.bedGraph $_output
rm -f ${n}_genome_merged.bedGraph

