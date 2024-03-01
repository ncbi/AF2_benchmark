## Overall steps
1. Download PDBs corresponding to both conformations of fold-switching proteins
2. Any proteins where one structure included only a short fragment or that had long gaps in the fold-switching region were excluded from the AF2Rank protocol. The final dataset consisted of 76 proteins (PDB IDs highlighted in supporting data of Table S1A)
3. To ensure that we passed the same sequence to AF2 for fold-switched conformations, we truncated extraneous N- and C-terminal residues used for purification but endogenous to their respective sequences.
4. If one structure included a domain that was not present in the other structure, that protein was excluded from the dataset
5. Any short gaps in the structures were modeled with RosettaCM, and the top scoring Rosetta model (minimum 1000 models generated) with a TM-score greater than 0.9 compared to the native structure were then selected for use
6. Hetero-atoms from non-standard residues such as the selenium in seleno-methionine and seleno-cysteine were replaced with their standard analogs (e.g. methionine and cysteine) using RosettaCM
7. Structures corresponding to each fold-switched conformation were passed to the AF2Rank protocol as templates, and the candidate structure’s accuracy is assessed based on confidence scores of the AF2 output model

## Important files/directories
1. gs_es.xlsx: table with all confidence scores generated from AF2Rank runs
2. runAF2Rankv4.sh: bash script to run AF2Rank protocol
3. prepare_pdb_chains.py: python script to process PDBs before being passed to AF2Rank
4. EAC_AF2Rank/: Fork of AF2Rank directory (https://github.com/jproney/AF2Rank) debugged to run on cluster
5. figures/: Figures of AF2 confidence scores. These include data for all fold-switching proteins in the non-excluded dataset
6. figures_exclusions/: Figures of AF2 confidence scores from the excluded dataset (Steps 2 and 4 of Overall steps, see above)
7. all_folds2_output/: Output directory generated from AF2Rank

For more information on running AF2Rank, see https://github.com/jproney/AF2Rank

References:
1. Roney, J. P., & Ovchinnikov, S. (2022). State-of-the-art estimation of protein model accuracy using AlphaFold. Physical Review Letters, 129(23), 238101.
2. Chakravarty, D., Schafer, J. W., Chen, E. A., Thole, J., & Porter, L. (2023). AlphaFold2 has more to learn about protein energy landscapes. bioRxiv, 2023-12.

