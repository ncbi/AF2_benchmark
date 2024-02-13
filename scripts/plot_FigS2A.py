import matplotlib.pyplot as plt
import scipy.stats as stats
from pathlib import Path
import numpy as np


def plot_PDF(x):
    # set plot parameters
    plt.rcParams.update({'font.size': 20})
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.family"] = "Helvetica"
    plt.figure(figsize=(12, 8))

    fit = stats.norm.pdf(x, np.mean(x), np.std(x))
    plt.plot(y, fit, '-', color='black', linewidth=2,
             alpha=0.75)
    plt.axvline(5, color='red', linestyle='--', linewidth=2)
    plt.xlabel('RMSD (Å) of Fold-switching region')
    plt.ylabel('Normalized frequency')
    plt.title("Data from AF-cluster runs w/o templates")
    plt.savefig(str(output_file), dpi=300)

# setting path and directories
main_dir = Path("/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo")
plot_dir = main_dir / "plots"
output_file = plot_dir / "FigS2A.png"
data_dir = main_dir / "data"

inputfile = data_dir / "rmsd_fs_all.npy"
y = np.load(inputfile)
y.sort()
vals2 = (y > 0.6).sum()
med2 = np.median(y)
print(vals2, y.size, round(np.average(y), 2), med2)
print("Top10 ", round((vals2/y.size), 2)*100)
plot_PDF(y)
