"""
plotting Heatmap
of TMscores for Sa1
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# setting figure parameters
plt.figure(figsize=(10, 14))
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Helvetica"

# setting path and directories
main_dir = Path("/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo")
plot_dir = main_dir / "plots"
output_file = plot_dir / "Fig3A.png"
data_dir = main_dir / "data"

# setting all background to transparent
plt.rcParams.update({
    "figure.facecolor":  (0.0, 0.0, 0.0, 0.0),  # red   with alpha = 30%
    "axes.facecolor":    (0.0, 0.0, 0.0, 0.0),  # green with alpha = 50%
    "savefig.facecolor": (0.0, 0.0, 0.0, 0.0),  # blue  with alpha = 20%
})

input_file = data_dir / "8e6y_sa1_v90t_30c/tmscores.txt"  # TM-scores file

# define the data list
out_data1 = {}
out_data2 = {}

with open(input_file, 'r') as oF:
    for line in oF:
        line = line.strip().strip(',')
        score_lst = line.split(':')
        run = score_lst[0].strip()
        scores = [tuple(float(c) for c in a.split())
                  for a in score_lst[1].split(',')]
        for score in scores:
            score1 = float(score[0])
            score2 = float(score[1])
            if run not in out_data1:
                out_data1[run] = [score1]
            else:
                out_data1[run].append(score1)
            if run not in out_data2:
                out_data2[run] = [score2]
            else:
                out_data2[run].append(score2)
my_ranks = []
for i in range(50):
    rank = "M" + str(i)
    my_ranks.append(rank)

# create the pandas dataframe
df1 = pd.DataFrame.from_dict(out_data1, orient='index', columns=my_ranks)
df2 = pd.DataFrame.from_dict(out_data2, orient='index', columns=my_ranks)


fig, ax = plt.subplots(1, 2, figsize=(12, 6), dpi=150)

sns.set(font_scale=1.5)

ax1 = sns.heatmap(df1, cmap='PiYG', yticklabels=False,
                  vmin=0, vmax=0.8, ax=ax[0], cbar=True)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=60, horizontalalignment='right',
                    fontweight='light',
                    fontsize='x-small')
ax2 = sns.heatmap(df2, cmap='PiYG', yticklabels=False,
                  vmin=0, vmax=0.8, ax=ax[1], cbar=True)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=60, horizontalalignment='right',
                    fontweight='light',
                    fontsize='x-small')
ax1.set_ylabel('TM-scores', fontsize = 20, fontweight='bold')
ax2.set_ylabel('TM-scores', fontsize = 20, fontweight='bold')
ax1.set_title("α/β-plait", fontsize=20, fontweight="bold")
ax2.set_title("3-helix bundle", fontsize=20, fontweight="bold")
ax1.invert_yaxis()
ax2.invert_yaxis()

plt.savefig(str(output_file), edgecolor='black', dpi=600, transparent=True)