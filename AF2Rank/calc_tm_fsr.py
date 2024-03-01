# calculate fold-switching region TM score and PLDDT
# 10-6-2023
import os
import sys
import pymol
from tmalign import tmscore
import pandas as pd
from Bio.PDB import PDBParser
from Bio import AlignIO
from Bio.SeqUtils import seq1
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore',BiopythonDeprecationWarning)
from Bio import pairwise2


TM_exe = "/Users/chenea/work/VSCode/aces/TMscore"

fsr_indices_csv = "/Users/chenea/work/range_fs_pairs_all.csv"
fsr_seq_xlsx = "/Users/chenea/work/VSCode/af2rank_0923/TableS1A.xlsx"
ground_state_xlsx = "/Users/chenea/work/pro4353-sup-0006-tables5.xlsx"

models_dir = "/Users/chenea/work/VSCode/af2rank_results/all_folds2_output/pdbs/"
results_dir = "/Users/chenea/work/VSCode/af2rank_results/all_folds2_output/results"
inputs_dir = "/Users/chenea/work/VSCode/af2rank_0923/final2_and_debug2_inputs"

def calc_tm(pdbid, model, idx):
    # calculate tm of fold switching region of a pdb
    #print("running for {:}".format(pdbid))
    model1 = models_dir+"{:}_decoy_{:}.pdb".format(pdbid,pdbid)
    template1 = "{:}/{:}.pdb".format(inputs_dir, pdbid)
    #print(model1)
    #print(template1)

    pymol.cmd.load(model1)
    pymol.cmd.load(template1)
    names = pymol.cmd.get_names()

    #score1 = tmscore('{:}'.format(names[0]),'{:}'.format(names[1]),quiet=0,exe=TM_exe)
    score1 = tmscore('{:} and i. {:}-{:}'.format(names[0],idx[0],idx[1]),'{:} and i. {:}-{:}'.format(names[1],idx[0],idx[1]),quiet=1,exe=TM_exe)
    pymol.cmd.delete('all')

    # calculate plddt
    bfactors = []
    for chain in model:
        for residue in chain:
            resid = residue.get_id()[1]
            #sys.exit()
            for atom in residue:
                if atom.name == "CA" and resid >= idx[0] and resid <= idx[1]:
                    bfactors.append(atom.bfactor)
    try:
        plddt = sum(bfactors)/len(bfactors)
    except:
        print("bfactors has length 0")
        sys.exit()

    return score1, plddt

def read_pdb(pdbid, pdb_path, parser):
    # read in the sequence of a pdb
    with warnings.catch_warnings(record=True) as w:  
        warnings.simplefilter("always", PDBConstructionWarning)
        structure = parser.get_structure(pdbid,pdb_path)
    model = structure[0]
    aa_seq = ""
    for chain in model:
        for residue in chain:
            aa_seq += seq1(residue.get_resname())
    return aa_seq, model

def align_fsr(pdb_seq,fsr_seq):
    a = pairwise2.align.localxs(pdb_seq,fsr_seq,-1,-0.5)
    i_start = None
    i_end = None
    leading_gaps = 0 # number of leading gaps of pdb_seq in alignment (correcting factor)
    read_status = {0: False, 1: False} 

    # find i_start
    for i in range(len(a[0][0])):

        # read statuses
        for j in read_status.keys():
            if not read_status[j] and a[0][j][i] != "-":
                read_status[j] = True
                if j == 0:
                    leading_gaps = i
        
        if read_status[0] and read_status[1]:
            i_start = i+1-leading_gaps
            break

    # find i_end
    read_status = {0: False, 1: False} 
    for i in reversed(range(len(a[0][0]))):

        # read statuses
        for j in read_status.keys():
            if not read_status[j] and a[0][j][i] != "-":
                read_status[j] = True
        
        if read_status[0] and read_status[1]:
            i_end = i+1-leading_gaps
            break

    return [i_start,i_end]

def read_results_excel(excel_path):
    # going to not use pandas, seems like its slower than just reading a file
    with open(excel_path, "r") as f:
        lines = f.readlines()
    decoy_info = lines[-1][:-1].split(",")
    #print(decoy_info)
    tm = eval(decoy_info[10])
    plddt = eval(decoy_info[12])
    ptm = eval(decoy_info[13])
    return tm, plddt, ptm

def get_ground_and_excited_states():

    ground_state_df = pd.read_excel(ground_state_xlsx)
    ground_states = {}

    for x in range(len(ground_state_df["Fold1"])):
        fold1 = ground_state_df["Fold1"][x]
        fold2 = ground_state_df["Fold2"][x]

        es = ""
        gs = ""

        # check which one is excited state
        if fold1[-2:] == "**" and fold2[-2:] != "**":
            es = fold1[:-2]
            gs = fold2
        elif fold2[-2:] == "**" and fold1[-2:] !="**":
            es = fold2[:-2]
            gs = fold1
        else:
            #print(fold1+" "+fold2)
            #continue
            #es = fold1.replace("*", "")
            #gs = fold2.replace("*", "")
            gs = "equilibrium"
        ground_states[fold1.replace("*","")+fold2.replace("*","")]=gs
    
    return ground_states

def main():

    print("reading excel files...")
    idx_df = pd.read_csv(fsr_indices_csv)
    fsr_seq_df = pd.read_excel(fsr_seq_xlsx)
    fsr_seq_dic = {}

    # get dictionary of fold_switchers and their fsrs
    for i in range(len(fsr_seq_df["Fold1"])):
        fold1_id = fsr_seq_df["Fold1"][i][:-1]+"_"+fsr_seq_df["Fold1"][i][-1:]
        fold2_id = fsr_seq_df["Fold2"][i][:-1]+"_"+fsr_seq_df["Fold2"][i][-1:]
        fsr_seq_dic[fold1_id] = fsr_seq_df["Sequence of fold-switching region"][i]
        fsr_seq_dic[fold2_id] = fsr_seq_df["Sequence of fold-switching region"][i]

    # get ground states
    ground_states_dic = get_ground_and_excited_states()
    ground_states = []

    fsr_tms = [[],[]]
    fsr_plddts = [[],[]]
    pdb_ids = [[], []]
    tms = [[], []]
    plddts = [[], []]
    ptms = [[], []]
    composites = [[],[]]

    print("parsing and scoring pdbs...")
    p = PDBParser()

    for i in range(len(idx_df["# pdb1"])):

        pdb1 = idx_df["# pdb1"][i]
        pdb2 = idx_df["pdb2"][i]
        
        model1_path = "{:}/{:}_decoy_{:}.pdb".format(models_dir, pdb1, pdb1)
        model2_path = "{:}/{:}_decoy_{:}.pdb".format(models_dir, pdb2, pdb2)
        if os.path.exists(model1_path) and os.path.exists(model2_path):            
            pdb1_seq, model1 = read_pdb(pdb1, model1_path, p) 
            pdb2_seq, model2 = read_pdb(pdb2, model2_path, p)

            #print(pdb1)
            idx1 = align_fsr(pdb1_seq,fsr_seq_dic[pdb1])
            idx2 = align_fsr(pdb2_seq,fsr_seq_dic[pdb2])


            # calculate tmscores and plddt scores
            fsr_tm1, fsr_plddt1 = calc_tm(pdb1, model1, idx1)
            fsr_tm2, fsr_plddt2 = calc_tm(pdb2, model2, idx2)

            tm1, plddt1, ptm1 = read_results_excel("{:}/results_{:}.csv".format(results_dir,pdb1))
            tm2, plddt2, ptm2 = read_results_excel("{:}/results_{:}.csv".format(results_dir,pdb2))

            fsr_tms[0].append(fsr_tm1)
            fsr_tms[1].append(fsr_tm2)
            fsr_plddts[0].append(fsr_plddt1)
            fsr_plddts[1].append(fsr_plddt2)
            pdb_ids[0].append(pdb1)
            pdb_ids[1].append(pdb2)
            tms[0].append(tm1)
            tms[1].append(tm2)
            plddts[0].append(plddt1)
            plddts[1].append(plddt2)
            ptms[0].append(ptm1)
            ptms[1].append(ptm2)
            composites[0].append(tm1*plddt1*ptm1)
            composites[1].append(tm2*plddt2*ptm2)
            try:
                ground_states.append(ground_states_dic[pdb1+pdb2])
            except:
                ground_states.append(ground_states_dic[pdb2+pdb1])
        #sys.exit()
    #tmscore('{:} and i. {:}-{:}'.format(names[0],fs_start,fs_stop),'{:} and i. {:}-{:}'.format(names[1],fs_start,fs_stop),quiet=0,exe=TM_exe)

    '''
    Still to do:
    DONE Calculate plddt for fold-switching region
    Package all metrics: tmscore, plddt, ptm, composite, tm(fsr), and ptm(fsr)
    into a single excel spreadsheet. From that I can make graphs with matplotlib
    Just going to test writing output to an excel sheet DONE

    All I need to do now is read in the individual excel sheets for each fold-switcher
    '''

    # write output
    excel_data = [pdb_ids[0],fsr_tms[0],fsr_plddts[0],tms[0],plddts[0],ptms[0],composites[0],
                  pdb_ids[1],fsr_tms[1],fsr_plddts[1],tms[1],plddts[1],ptms[1],composites[1], ground_states]
    excel_writer = list(map(list, zip(*excel_data)))
    output_df = pd.DataFrame(excel_writer,
                             columns=["pdb1", "fsr_tm1", "fsr_plddt1", "tm1", "plddt1", "ptm1", "composite1",
                                      "pdb2", "fsr_tm2", "fsr_plddt2", "tm2", "plddt2", "ptm2", "composite2", "ground_states"])
    output_df.to_excel("af2rank_output_data.xlsx")

    '''for i in range(len(fsr_tms[0])):
        print("{:} {:} {:} {:} {:} {:}".format(pdb_ids[0][i],pdb_ids[1][i], fsr_tms[0][i], fsr_tms[1][i], fsr_plddts[0][i], fsr_plddts[1][i]))'''

if __name__ == "__main__":
    main()
    #read_results_excel("/Users/chenea/work/VSCode/af2rank_results/all_folds2_output/results/results_4dxt_A.csv")
    #get_ground_and_excited_states()