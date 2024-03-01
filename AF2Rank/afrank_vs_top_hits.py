import pandas as pd
afrank_df = pd.read_excel("/Users/chenea/work/VSCode/af2rank_0923/af2rank_output_data.xlsx")
all_hits_df = pd.read_excel("/Users/chenea/work/VSCode/af2rank_0923/gs_es.xlsx", "all")
top_hits_df = pd.read_excel("/Users/chenea/work/VSCode/af2rank_0923/gs_es.xlsx", "high_quality")

fold1_all_hits = all_hits_df["Fold1All"]
fold2_all_hits = all_hits_df["Fold2All"]
fold1_high_hits = top_hits_df["Fold1All"]
fold2_high_hits = top_hits_df["Fold2All"]
pdb1s = top_hits_df["pdb1"]
pdb2s = top_hits_df["pdb2"]

for i in range(len(pdb1s)):
    if fold1_all_hits[i] != 0 and fold1_high_hits[i] == 0:
        try:
            j = list(afrank_df["pdb1"]).index(pdb1s[i])
        except:
            continue
        score1 = afrank_df["fsr_tm1"][j]
        score2 = afrank_df["fsr_tm2"][j]
        print("fold1 missing {:} {:} fsr_tm scores: {:} {:}".format(pdb1s[i],pdb2s[i],score1,score2))
    if fold2_all_hits[i] != 0 and fold2_high_hits[i] == 0:
        try:
            j = list(afrank_df["pdb2"]).index(pdb2s[i])
        except:
            continue
        score1 = afrank_df["fsr_tm1"][j]
        score2 = afrank_df["fsr_tm2"][j]
        print("fold2 missing {:} {:} fsr_tm scores: {:} {:}".format(pdb1s[i],pdb2s[i],score1,score2))