"""
plotting Heatmap
of TMscores for BCCIP variants
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def plot_heatmap(inputfile,title,outputfile):
    """
    read data
    create pandas dataframe
    plot and save heatmap
    """
    # setting figure parameters
    plt.figure(figsize=(5, 5))
    plt.rcParams.update({'font.size': 20})
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.family"] = "Helvetica"
    
    # setting all background to transparent
    plt.rcParams.update({
        "figure.facecolor":  (0.0, 0.0, 0.0, 0.0),  # red   with alpha = 30%
        "axes.facecolor":    (0.0, 0.0, 0.0, 0.0),  # green with alpha = 50%
        "savefig.facecolor": (0.0, 0.0, 0.0, 0.0),  # blue  with alpha = 20%
    })
    
    out_dict = {}
    with open(inputfile) as Infile:
        for line in Infile:
            line=line.strip('\n')
            line=line.strip(',')
            tmscores = []
            (key,rest)=line.split(':')
            val_lst = rest.split(',')
            tmscores = [float(x) for x in val_lst]
            if key not in out_dict:
                out_dict[key]=tmscores

    my_ranks = []
    for i in range(50):
        rank = "M" + str(i)
        my_ranks.append(rank) 

    df = pd.DataFrame.from_dict(out_dict, orient='index',columns=my_ranks)

    sns.set(font_scale=1.5)

    ax = sns.heatmap(df,cmap='PiYG',yticklabels=False,vmin=0, vmax=0.8)
    ax.invert_yaxis()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=60,horizontalalignment='right', 
                   fontweight='light',
                   fontsize='x-small' )

    plt.ylabel('TM-scores', fontsize = 20,fontweight="bold")
    ax.set_title(title, fontsize=20, fontweight="bold")

    plt.savefig(str(outputfile), edgecolor='black', dpi=600, transparent=True)
    plt.close()

# setting path and directories
main_dir = Path("/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo")
plot_dir = main_dir / "plots"
data_dir = main_dir / "data"

input_file1 = data_dir / "bccip_vars/tmscores_b.txt"  # TM-scores file for BCCIPα
input_file2 = data_dir / "bccip_vars/tmscores_a.txt"  # TM-scores files for BCCIPβ
output_file1 = plot_dir / "Fig3B_leftpanel.png"
output_file2 = plot_dir / "Fig3B_rightpanel.png"

plot_heatmap(input_file1, "BCCIPβ",output_file1)
plot_heatmap(input_file2, "BCCIPα",output_file2)

