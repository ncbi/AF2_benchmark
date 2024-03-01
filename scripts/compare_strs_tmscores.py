"""
Compare the predicted models with original PDBs
report TM-scores for ranked 0 to 4
input line is pdb1 pdb2 predicted_models_path

This version requires tmtools 0.0.2 (Python bindings around the TM-align code for structural alignment of proteins)
check this for local installation
https://pypi.org/project/tmtools/

Usage:

python3.8 compare_strs_tmscores.py pdb1 pdb2 <path_to_preds>

"""
import sys
from pathlib import Path
import glob
# call related modules of tmtools after installation
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path


def get_tmscore(r,predfilepath):
	"""
	parameters:
	r - tmtool obj for original PDB
	predfilepath - path for predicted files

	returns:
	tmscore list

	"""
	
	coords2, seq2 = get_residue_data(r)
	tmscores = []
	modelfiles = sorted(glob.glob(str(predfilepath) + "/*_relaxed*pdb"))

	for model in modelfiles:
		modelpath = Path(model)
		model  = str(modelpath.parent) + "/" + modelpath.stem 
		s = get_structure(get_pdb_path(model))
		coords1, seq1 = get_residue_data(s)
		res = tm_align(coords1, coords2, seq1, seq2)
		tmscore = round(res.tm_norm_chain1,2) # wrt to model
		tmscores.append(tmscore)

	return tmscores

def run_for_models(FH,data_dir,r1,r2):
	"""
	compare the original PDB
	with the predicted models ranked 0 to 4
	
	parameters:
	FH - filehandle for writing
	data_dir - path for the predicted strs
	r1,r2 - tmtool objects for pdb1 and pdb2
	
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

		tmscore_lst1 = get_tmscore(r1,preddir) # wrt pdb1
		tmscore_lst2 = get_tmscore(r2,preddir) # wrt pdb2
		
		for index in range(0,5):
			outline2 = str(tmscore_lst1[index]) + " " + str(tmscore_lst2[index]) + ","
			FH.write(outline2)
		FH.write("\n")
			

# setting path
data_dir = Path(
	'/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data')
pdb_dir = data_dir / "clean_chains"
pred_dir = data_dir 

# input arguments
pdb1 = sys.argv[1]
pdb2 = sys.argv[2]
inputpath = pred_dir / sys.argv[3]

print("Running for pair ", pdb1, pdb2, end="...")

pdbfile1 = pdb_dir / pdb1
pdbfile2 = pdb_dir / pdb2
# getting the coords for the original PDBs
r1 = get_structure(get_pdb_path(str(pdbfile1)))
r2 = get_structure(get_pdb_path(str(pdbfile2)))

# open files for saving all tmscores as output
# in requisite folders
outfilepath = Path (inputpath) /  "tmscores.txt"

OF = open(outfilepath,'w')

run_for_models(OF,inputpath,r1,r2) # predicted models compared with PDB1 and PDB2

OF.close()

print("Done")
