# there are some inputs in the debug2 folder, and others in the final2 folder.
# combine both folders. For pdbs that exist in both folders, use the one from the
# debug2 directory
import os
import sys

final2_dir = "/Users/chenea/work/RankProject/all_folds2/inputs/final2"
debug2_dir = "/Users/chenea/work/RankProject/all_folds2/inputs/debug2"
finished_targets_file = "/Users/chenea/work/RankProject/all_folds2/all_folds2_output/finished_targets.txt"
new_dir = "/Users/chenea/work/VSCode/af2rank_0923/final2_and_debug2_inputs"

def get_pdb_paths(pdb_dir, pdb_dic):
    new_dic = pdb_dic

    for f in os.listdir(pdb_dir):
        if f.endswith(".pdb"):
            new_dic[f[:-4]] = pdb_dir+"/"+f

    return new_dic

def main():

    # read in finished targets
    finished_targets = []
    with open(finished_targets_file, "r") as f:
        finished_targets = [line[:-1] for line in f]

    # read in debug dir
    pdb_dic = {}
    pdb_dic = get_pdb_paths(final2_dir, pdb_dic)
    pdb_dic = get_pdb_paths(debug2_dir, pdb_dic)

    # check that all finished targets are in pdb_dic, and vice versa
    for target in finished_targets:
        if target not in pdb_dic.keys():
            print("{:} not in pdb_dic.keys()".format(target))
    
    print("copying files to new directory: {:}".format(new_dir))
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    for key in pdb_dic.keys():
        if key not in finished_targets:
            print("{:} not in finished targets".format(key))
        else:
            os.system("cp {:} {:}/{:}.pdb".format(pdb_dic[key], new_dir, key))
            #sys.exit()

    
main()

