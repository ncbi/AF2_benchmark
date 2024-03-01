# Scripts for efficiently generating the sets of contacts to be mutated to Alanine by the SPEACH-AF algorithm


## To run the code:
## 1.  tar -xzf clean_pdbs.tgz
## 2.  python /get_contact_sets.py 1cee_test.csv clean_pdbs 1cee_contacts
##
### This will generate all combinations of residues that will be mutated to alanine in the MSA.  Each number corresponds to its respective column in the MSA, starting from 0.
### The result should match 1cee_contacts.tgz
###
### range_fs_pairs_all.txt contains all fold switch pairs and their fold-switching regions (ranges of residues sampled for contacts to be mutated to alanine).  The third and fourth columns are PDB numbers (not used but provided for reference); fifth and sixth are residue indices (used to get contacts)
### For more information about how and why these mutations are chosen, see: Stein and Mchaourab (2022) PLOS Comp. Bio.
