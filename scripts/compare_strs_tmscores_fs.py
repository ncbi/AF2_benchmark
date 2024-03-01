"""
Compare the predicted models with original PDBs
report TM-scores for ranked 0 to 4
input line is pdb1 pdb2 preds_of_pdb dirname

Requires tmtools 0.0.2 (Python bindings around the TM-align code for structural alignment of proteins)
check this for local installation
https://pypi.org/project/tmtools/

Usage:

python3.8 compare_strs_tmscores_fs.py pdb1 pdb2 preds <path_to_preds>

"""
import sys
from pathlib import Path
import glob

# call related modules of tmtools after installation
from tmtools import tm_align
import numpy as np
from Bio.PDB import PDBParser

pdbParser = PDBParser(QUIET=True)

# convert three letter code to one letter code
aa3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# get numpy arrays for coords at the fold-switching region
# also return the seq in 1-letter code for the same

def get_coords(pdbfile,fs_range):
        """
        parameters:
        pdbfile - path to pdbfile
        fs_range - range of residues at the fold-switching region, given as string - "112-162"
        returns:
        numpy array of coords
        string of seqs in 1-letter-code
        """

        seq = ""
        struct = pdbParser.get_structure('x',str(pdbfile))
        coords = []
        seq_dict = {}

        # for residues within a certain range, using numpy to save the coords
        # and save the sequence as a dict and then sorted list of tuples
        # return the coords and the seq

        # convert str to residue range for the fs region
        (start,stop) = fs_range.split("-")
        res_range = range(int(start),int(stop)+1)

        for atom in struct.get_atoms():
                residue = atom.get_parent() # from atom we can get the parent residue
                res_id = residue.get_id()[1]
                resname = residue.get_resname()
                if res_id in res_range and atom.get_name()=="CA":
                        x,y,z = atom.get_coord()
                        coords.append([x,y,z])
                        if res_id not in seq_dict:
                                seq_dict[res_id]=aa3to1[resname]
 
        # convert to np array
        coords_np = np.array(coords)
        # sort the seq_dict by keys a.k.a res_ids
        sorted_data = sorted(seq_dict.items())
        for i in sorted_data:
                seq+=i[1]

        return  coords_np,seq


def get_tmscore(coords1,seq1,predfilepath,res_range):
	"""
	parameters:
	coords1, seq1 - the numpy array of PDB coords and its seqs
	predfilepath - path for predicted files
	res_range - fs range in predicted models

	returns:
	tmscore list

	"""
	
	tmscores = []
	modelfiles = sorted(glob.glob(str(predfilepath) + "/*_relaxed*pdb"))

	for model in modelfiles:
		modelpath = Path(model)
		coords2, seq2 = get_coords(modelpath,res_range)
		res = tm_align(coords1, coords2, seq1, seq2)
		tmscore = round(res.tm_norm_chain1,2) # wrt to model
		tmscores.append(tmscore)

	return tmscores

def run_for_models(FH,pdbfile1,pdbfile2,data_dir,pred_range,res_range1,res_range2):
	"""
	compare the original PDB
	with the predicted models, 0 to 5

	parameters:
	FH - filehandle for writing
	pdbfile1 - path to original PDB, Fold1
	pdbfile2 - path to alternate PDB, Fold2
	data_dir - path for the predicted strs
	res_range1 - fs range in PDB1 and its models
	res_range2 - fs range in PDB2 and its models
        
	returns:
	nothing
        
	saves the TM-scores in a local file
	"""
	# get list of subdirectories
	all_sub_dir_paths = glob.glob(str(data_dir) + '/*/') # returns list of sub directory paths
	# files found then continue	
	if len(all_sub_dir_paths) == 0:
		pass
	for subdir in all_sub_dir_paths:
		filename = Path(subdir).name
		outline1 = filename + " : "
		FH.write(outline1)
		
		preddir = Path(subdir)

		# predicted dir doesn't exist then continue
		if not preddir.exists():
			pass
		# only comparing on one set of predicted models
        	# but with both PDBs/Folds
		coords1,seq1 = get_coords(pdbfile1,res_range1)
		coords2,seq2 = get_coords(pdbfile2,res_range2)
		tmscore_lst1 = get_tmscore(coords1,seq1,preddir,pred_range) # wrt pdb1
		tmscore_lst2 = get_tmscore(coords2,seq2,preddir,pred_range) # wrt pdb2
	
		for index in range(0,5):
			outline2 = str(tmscore_lst1[index]) + " " + str(tmscore_lst2[index]) + ","
			FH.write(outline2)
		FH.write("\n")
			

# setting path
data_dir = Path(
    '/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data')
pdb_dir = data_dir / "clean_chains"
pred_dir = data_dir 

# the range of the fold-switching region
range_file = data_dir / "range_fs_pairs_test.txt"

# convert this file into a dictionary for reference later
fs_res = {}

# The range_file file has the fold-switching residue ranges
# for the original PDB/PDB1, alternate PDB/PDB2
# Predicted model for PDB1, predicted model for PDB2
with open(range_file, 'r') as Infile:
    next(Infile)  # skip header line "# pdb1,pdb2,pred1,pred2"
    for line in Infile:
        line = line.strip()
        (n1, n2, p1, p2, m1, m2) = line.split(",")
        # the value of the dictionary is a tuple
        # the first element of tuple is the fs range in the original PDB
        # followed by the range in the predicted model
        if n1 not in fs_res:
            fs_res[n1] = (p1, m1)
        if n2 not in fs_res:
            fs_res[n2] = (p2, m2)

# input arguments
pdb1 = sys.argv[1]  # 6c6s_D
pdb2 = sys.argv[2]  # 2oug_C
preds = sys.argv[3]  # 2oug_C
inputpath = sys.argv[4]  # 2oug_C

data_dir = pred_dir / inputpath  # Path to the predicted models

print("Running for pair ", pdb1, pdb2, end="...")
print("evaluating predictions of ", preds, end="...")

# getting residue ranges of fold-switching regions
# say we are running for pair 1nqd_A,1nqj_B
try:
    # so if pdb1 is '1nqd_A', fs_res['1nqd_A']=('895-919', '1-33')
    range_pdb1 = fs_res[pdb1]
    # and if pdb2 is '1nqj_B', fs_res['1nqj_B']=('894-919', '1-33')
    range_pdb2 = fs_res[pdb2]
except:
    print("check PDBIDs ", pdb1, pdb2)
    sys.exit(1)

if preds == pdb1:
    range_pred = range_pdb1[1]
else:
    range_pred = range_pdb2[1]

pdbfile1 = pdb_dir / f"{pdb1}.pdb"
pdbfile2 = pdb_dir / f"{pdb2}.pdb"


# open files for saving all tmscores as output
# in requisite folders
outfilepath = data_dir / "tmscores_fs.txt"

OF = open(outfilepath, 'w')

# predicted models compared with PDB1 and PDB2
run_for_models(OF, pdbfile1, pdbfile2, data_dir, range_pred, range_pdb1[0], range_pdb2[0])

OF.close()

print("Done")
