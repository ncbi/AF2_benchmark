"""
Requires align module written locally

Usage:
python3.8 compare_strs_fs.py pdb1 pdb2 preds <path_to_preds>
"""

import sys
from align import cal_rmsd_partial
from pathlib import Path
import glob


def get_rmsd(predfilepath,pdbfilepath,start1,stop1,start2,stop2):
        """
        parameters:
        pdbfilepath - path for original PDB file
        predfilepath - path for predicted files
        
        start,stop indices for fs region

        returns:
        rmsd list

        """
        rmsd_vals = []
        modelfiles = sorted(glob.glob(str(predfilepath) + "/*_relaxed*pdb")) 
        
        for model in modelfiles:
                modelpath = Path(model)
                try:
                    rmsd = round(cal_rmsd_partial(modelpath,pdbfilepath,start1,stop1,start2,stop2),2)
                    rmsd_vals.append(rmsd)
                except AssertionError as error:
                        print(f"NOTE: CANNOT align {modelpath} {pdbfilepath} because of {error}")

        return rmsd_vals

def run_for_models(FH,pdbfile1,pdbfile2,preddir,pred_range,res_range1,res_range2):
        """
        compare the original PDB
        with the predicted models ranked 0 to 5

        parameters:
        FH - filehandle for writing
        preddir - path for the predicted strs
        pdbfile1,pdbfile2 - full paths for pdb1 and pdb2

        res_range1 - fs range in PDB1 and its models
        res_range2 - fs range in PDB2 and its models
        
        returns:
        nothing

        saves the RMSD values in a local file

        """

        # get start, stop indices for the fs region
        (start,stop) = pred_range.split("-") # model
        (start1,stop1) = res_range1.split("-") # pdb1
        (start2,stop2) = res_range2.split("-") # pdb2
        
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
                rmsd_lst1 = get_rmsd(preddir,pdbfile1,start,stop,start1,stop1) # wrt pdb1
                rmsd_lst2 = get_rmsd(preddir,pdbfile2,start,stop,start2,stop2) # wrt pdb2

                for index in range(0,5):
                        outline2 = str(rmsd_lst1[index]) + " " + str(rmsd_lst2[index]) + ","
                        FH.write(outline2)
                FH.write("\n")


# setting path
data_dir = Path('/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data')
pdb_dir = data_dir / "clean_chains"
pred_dir = data_dir 

# the range of the fold-switching region
range_file = data_dir / "range_fs_pairs_test.txt"

# convert this file into a dictionary for reference later
fs_res = {}

# The range_file file has the fold-switching residue ranges
# for the original PDB/PDB1, alternate PDB/PDB2
# Predicted model for PDB1, predicted model for PDB2
with open(range_file,'r') as Infile:
        next(Infile) # skip header line "# pdb1,pdb2,pred1,pred2"
        for line in Infile:
                line=line.strip()
                (n1,n2,p1,p2,m1,m2)=line.split(",")
                # the value of the dictionary is a tuple
                # the first element of tuple is the fs range in the original PDB 
                # followed by the range in the predicted model
                if n1 not in fs_res:
                        fs_res[n1]=(p1,m1)
                if n2 not in fs_res:
                        fs_res[n2]=(p2,m2)

# input arguments
pdb1 = sys.argv[1] # 6c6s_D
pdb2 = sys.argv[2] # 2oug_C
preds = sys.argv[3] # 2oug_C
inputpath = sys.argv[4] # 2oug_C

data_dir = pred_dir / inputpath # Path to the predicted models

print("Running for pair ",pdb1,pdb2,end="...")
print("evaluating predictions of ", preds,end="...")

# getting residue ranges of fold-switching regions
# say we are running for pair 1nqd_A,1nqj_B
try:
        range_pdb1 = fs_res[pdb1] # so if pdb1 is '1nqd_A', fs_res['1nqd_A']=('895-919', '1-33')
        range_pdb2 = fs_res[pdb2] # and if pdb2 is '1nqj_B', fs_res['1nqj_B']=('894-919', '1-33')
except:
        print("check PDBIDs ",pdb1,pdb2)
        sys.exit(1)

if preds==pdb1:
        range_pred = range_pdb1[1]
else:
        range_pred = range_pdb2[1]

pdbfile1 = pdb_dir / f"{pdb1}.pdb"
pdbfile2 = pdb_dir / f"{pdb2}.pdb"


# open files for saving all rmsd as output
# in requisite folders
outfilepath = data_dir /  "rmsds_fs.txt"

OF = open(outfilepath,'w')

run_for_models(OF,pdbfile1,pdbfile2,data_dir,range_pred,range_pdb1[0],range_pdb2[0]) # predicted models compared with PDB1 and PDB2

OF.close()

print("Done")
