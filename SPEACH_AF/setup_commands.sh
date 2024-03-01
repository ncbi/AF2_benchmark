#This script was used to mutate MSAs to alanine using .pik files
#using SPEACH_AF_contacts/get_contact_sets.py

#! /bin/bash

for i in `ls *pik`; do test=$(cut -d'.' -f1 <<< $i);
		       mkdir $test; mkdir ${test}/msas;
		       cp /data/Chakravarty_llp/AlphaFold_notemplates/all_FS_notemplates/${1}*.fasta ${test}.fasta; 
		       cp /data/Chakravarty_llp/AlphaFold_notemplates/all_FS_notemplates/${1}*/features.pkl ${test} ; 
		       cp /data/Chakravarty_llp/AlphaFold_notemplates/all_FS_notemplates/${1}*/msas/pdb_hits.hhr ${test}/msas; 
		       ./mutate_stockholm.py /data/Chakravarty_llp/AlphaFold_notemplates/all_FS_notemplates/${1}*/msas/uniref90_hits.sto $i > ${test}/msas/uniref90_hits.sto; 
		       ./mutate_stockholm.py /data/Chakravarty_llp/AlphaFold_notemplates/all_FS_notemplates/${1}*/msas/mgnify_hits.sto $i > ${test}/msas/mgnify_hits.sto; 
		       ./mutate_a3m.py /data/Chakravarty_llp/AlphaFold_notemplates/all_FS_notemplates/${1}*/msas/bfd_uniclust_hits.a3m $i > ${test}/msas/bfd_uniclust_hits.a3m ; done

ls *fasta > fastas.txt;
../../AF2_run.py fastas.txt $PWD;
chmod u+x *af2.sh
#for i in `ls *af2.sh`; do sbatch --cpus-per-task=6 --partition=gpu --mem=40g --gres=gpu:v100x:1,lscratch:100 --time=04:00:00 $i; done
