"""
Requires align module written locally

Usage:
python3.8 compare_strs_rmsd.py pdb1 pdb2 <path_to_preds>
"""

import sys
from align import cal_rmsd
from pathlib import Path
import glob


def get_rmsd(predfilepath,pdbfilepath):
        """
        parameters:
        pdbfilepath - path for original PDB file
        predfilepath - path for predicted files

        returns:
        rmsd list

        """
        rmsd_vals = []
        modelfiles = sorted(glob.glob(str(predfilepath) + "/*_relaxed*pdb"))
        
        for model in modelfiles:
                modelpath = Path(model)
                try:
                    rmsd = round(cal_rmsd(modelpath,pdbfilepath),2)
                    rmsd_vals.append(rmsd)
                except AssertionError as error:
                        print(f"NOTE: CANNOT align {modelpath} {pdbfilepath} because of {error}")
                        rmsd = 20.0
                        rmsd_vals.append(rmsd)

        return rmsd_vals

def run_for_models(FH,data_dir,pdbfilepath1,pdbfilepath2):
        """
        compare the original PDB
        with the predicted models ranked 0 to 5

        parameters:
        FH - filehandle for writing
        preddir - path for the predicted strs
        pdbfilepath1,pdbfilepath2 - full paths for pdb1 and pdb2

        returns:
        nothing

        saves the RMSD values in a local file

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

                rmsd_lst1 = get_rmsd(preddir,pdbfilepath1) # wrt pdb1
                rmsd_lst2 = get_rmsd(preddir,pdbfilepath2) # wrt pdb2

                for index in range(0,5):
                        outline2 = str(rmsd_lst1[index]) + " " + str(rmsd_lst2[index]) + ","
                        FH.write(outline2)
                FH.write("\n")


# setting path
data_dir = Path('/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo/data')
pdb_dir = data_dir / "clean_chains"
pred_dir = data_dir

# input arguments
pdb1 = sys.argv[1]
pdb2 = sys.argv[2]
inputpath = pred_dir / sys.argv[3]

print("Running for pair ",pdb1,pdb2,end="...")

pdbfile1 = pdb_dir / f"{pdb1}.pdb"
pdbfile2 = pdb_dir / f"{pdb2}.pdb"

# open file for saving all rmsd values as output

outfilepath = Path (inputpath) /  "rmsds.txt"

OF = open(outfilepath,'w')

run_for_models(OF,inputpath,pdbfile1,pdbfile2) # predicted models of PDB1, compared with PDB1 and PDB2

OF.close()

print("Done")
