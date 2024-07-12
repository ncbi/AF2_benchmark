import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
 
# setting path and directories
main_dir = Path("/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo")
plot_dir = main_dir / "plots"
output_file = plot_dir / "Fig1A_revised.png"

# setting figure parameters
plt.figure(figsize=(10, 14))
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Helvetica"

# values are taken from success_with_TMscore_metric_revised.xlsx
val1 = [2, 9, 5, 2, 1, 7, 12, 36] # unique
val2 = [6, 2, 7, 5, 6, 11, 20, 20] # common

names = ["AF2.2.0", "AF2.3.1","AF2_multimer","AF3","SPEACH_AF","AF-cluster","All_AF","ACE"]

totals = [x + y for x, y in zip(val1, val2)]

df = pd.DataFrame({'unique':val1,'common':val2},index=names)
threshold = 92
ax = df.plot(kind='bar', stacked=True, color=['grey', 'black'])
ax.set_ylim(1, 100)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45,
                   fontweight='bold', horizontalalignment='right')
for i, v in enumerate(totals):
   ax.text(i, v + 0.2, str(v), fontsize=16, weight='bold', ha='center')
plt.axhline(threshold, color='red', ls='dashed', linewidth=2.0)
plt.ylabel("#Successes", fontsize=16, fontweight="bold")
legend = plt.legend(loc="best", edgecolor=(
    0, 0, 0, 1.), facecolor=(1, 1, 1, 0.0))
legend.get_frame().set_alpha(None)

# for TM-score metric
plt.savefig(output_file, bbox_inches='tight',
 edgecolor='black', dpi=1200, transparent=True)
