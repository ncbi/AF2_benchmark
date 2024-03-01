## Overall steps
1.	Create a file of FASTA seqs of interest. 
2.	Generate MSAs using ColabFold 1.5.2 (https://hpc.nih.gov/apps/colabfold.html) 
3.	Run AF-cluster to generate clusters and shallow MSAs.
4.	Feed the shallow MSAs to ColabFold generate predictions.

## example Bash commands/scripts for the steps
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
