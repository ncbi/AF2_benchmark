# Scripts for efficiently generating the sets of Contacts to be mutated by Alanine by the SPEACH-AF algorithm


## To run the code:
## 1.  tar -xzf clean_pdbs.tgz
## 2.  python /get_contact_sets.py 1cee_test.csv clean_pdbs 1cee_contacts
##
### This will generate all combinations of residues that will be mutated to alanine in the MSA.  Each number corresponds to its respective column in the MSA, starting from 0.
### The result should match 1cee_contacts.tgz
### For more information about how and why these mutations are chosen, see: Stein and Mchaourab (2022) PLOS Comp. Bio.
