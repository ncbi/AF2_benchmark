import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# setting path and directories
main_dir = Path("/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo")
plot_dir = main_dir / "plots"
output_file = plot_dir / "FigS2B.png"

plt.figure(figsize=(10,14))
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Helvetica"

# the values are taken from data/success_with_RMSD_metric.xlsx
values = [8, 10, 7, 6, 14, 25, 56]

names = ["AF2.2.0", "AF2.3.1","AF2_multimer","SPEACH_AF","AF-cluster","All_AF2","ACE"]
threshold=93

ax = sns.barplot(x=names, y=values,color='black')
ax.set_ylim(1, 100)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontweight='bold', horizontalalignment='right' )
for i, v in enumerate(values):
   ax.text(i, v + 0.2, str(v), fontsize=21, weight='bold', ha='center')

plt.axhline(threshold, color='red', ls='dashed',linewidth=2.0)
plt.ylabel("#SUCCESS",fontsize=20, fontweight="bold")
plt.savefig(str(output_file), edgecolor='black', dpi=1200, transparent=True)