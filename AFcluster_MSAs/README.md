## Overall steps
1.	Create a file of FASTA seqs of interest. 
2.	Generate MSAs using ColabFold 1.5.2 (https://hpc.nih.gov/apps/colabfold.html) 
3.	Run AF-cluster to generate clusters and shallow MSAs.
4.	Feed the shallow MSAs to ColabFold generate predictions.

## Example Bash commands/scripts for the steps
Step1 –
Create a batch input file (e.g., colabfold-search.sh) to create the MSAs:
"#!/bin/bash
module load colabfold
colabfold_search --threads $SLURM_CPUS_PER_TASK \
    all_vars.fa $COLABFOLD_DB rfah_vars"

Run the batch file.
sbatch --cpus-per-task=16 --mem=128g --gres=lscratch:100 colabfold-search.sh

Step2 –
Run AF-cluster on the MSAs generated-

python AF_Cluster/scripts/ClusterMSA.py 0 -i 0.a3m -o 0_msas

Step3-
Run Colabfold on the shallow MSAs generated from AF-cluster-
Create a batch input file to run ColabFold

#!/bin/bash

module load colabfold
colabfold_batch --num-relax 50 --num-seeds 10 --amber --use-gpu-relax 0_msas 0_msas_models/

Run the batch file.
sbatch --time=48:00:00 --cpus-per-task=8 --mem=48g --gres=lscratch:100,gpu:a100:1 --partition=gpu run_colabfold_0.sh

## Successful predictions generated from AF_cluster added 
Please combine and uncompress the these files as follows - 

"cat afcluster_success.tar.gz.aa afcluster_success.tar.gz.ab afcluster_success.tar.gz.ac afcluster_success.tar.gz.ad afcluster_success.tar.gz.ae | tar xvfz -"

The resulting folder will have details on the protein pairs as a .CSV file and sub-folders with PDB structures (predictions), PyMol sessions of the successful predictions.

## TM-scores and prediction confidence scores added for successful predictions
The compressed file afcluster_success_tmscores.tgz has the same folder format as the one having successful predictions (mentioned above). 

The predictions are ranked by confidence scores (percentage of residues with pLDDT > 70) followed by TM-score calculated using the fold-switching region. 

To generate the predictions with high confidence and good TM-scores (> 0.6), please use the MSAs provided in AFcluster_MSAs. Please uncompress the file afcluster_success_tmscores.tgz to get the folder "success_with_AFcluster_scores/"

For, e.g., to generate the best-ranked prediction for AIAT, 3t1p, see the first line in the AIAT/3t1p_A_tmscores_fs_all.csv

Description of header # name of file rank_pred, confidence score, TM-score1 (w.r.t Fold1), TM-score2 (w.r.t Fold2)
"7_047 001,98.11,0.4,0.67"

To get the MSAs of the best prediction, please use the sub_5/7_msas/ in AFcluster_MSAs (check info_all_runs,txt for that information) and use the 7_047.a3m file in that folder, submit that MSA file to ColabFold to generate predictions. 
