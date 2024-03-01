#!/bin/sh

source myconda
mamba activate AF2Rank_python-3o7
module load cuDNN/7.6.5/CUDA-10.1

WKDIR="/data/Chen_SIP_2022/PostBac23/RankProject/EAC_AF2Rank"
name="all_folds2_output"
target="/data/Chen_SIP_2022/PostBac23/RankProject/all_folds2/all_foldsL_targets.txt"
#target_list="4n9w_A"
seed="1"
rec="1"
decoy_dir="/data/Chen_SIP_2022/PostBac23/RankProject/all_folds2/all_folds2_decoys/"
af2_dir="/data/Chen_SIP_2022/PostBac23/RankProject/alphafold/"
tm_exec="/data/Chen_SIP_2022/PostBac23/RankProject/TMscore"
output_dir="/data/Chen_SIP_2022/PostBac23/RankProject/all_folds2"

# Make directory
mkdir $output_dir/$name
mkdir $output_dir/$name/pdbs
mkdir $output_dir/$name/results

output_dir="/data/Chen_SIP_2022/PostBac23/RankProject/all_folds2/"

python $WKDIR/test_templates_debug1o5.py $name --targets_file $target --seed $seed --recycles $rec --decoy_dir $decoy_dir --seq_replacement - --mask_sidechains_add_cb --af2_dir $af2_dir --tm_exec $tm_exec --output_dir $output_dir
