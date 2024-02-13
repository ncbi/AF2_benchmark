import matplotlib.pyplot as plt
import scipy.stats as stats
from pathlib import Path
import numpy as np


def plot_PDF(x, y,title,outputfile):
    # set plot parameters
    plt.rcParams.update({'font.size': 20})
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.family"] = "Helvetica"
    plt.figure(figsize=(12, 8))
    fit1 = stats.norm.pdf(x, np.mean(x), np.std(x))
    fit2 = stats.norm.pdf(y, np.mean(y), np.std(y))
    plt.plot(x, fit1, '--', color='grey',
             linewidth=2, alpha=0.75, label='Whole')
    plt.plot(y, fit2, '-', color='black', linewidth=2,
             alpha=0.75, label='Fold-switching region')
    plt.xlabel('TM-scores')
    plt.ylabel('Normalized frequency')
    plt.title(title)
    plt.legend(loc='best')
    plt.savefig(str(outputfile), dpi=300)
    plt.close()

# setting path and directories
main_dir = Path(
    "/Users/chakravartyd2/fold_switching/mchaourab_project/for_github_repo")
plot_dir = main_dir / "plots"
output_file1 = plot_dir / "FigS1_distribution1.png"
output_file2 = plot_dir / "FigS1_distribution2.png"
data_dir = main_dir / "data"


# plot for AF2.3.1 data
inputfile1 = data_dir / "all_tmscores_af231_not.npy"
inputfile2 = data_dir / "all_fstmscores_af231_not.npy"
x = np.load(inputfile1)
x.sort()
y = np.load(inputfile2)
y.sort()
plot_PDF(x, y, 'Data from AF2.3.1 w/o templates', output_file1)


# plot for AF-cluster data
inputfile1 = data_dir / "tmscores_all.npy"
inputfile2 = data_dir / "tmscores_fs_all.npy"
x = np.load(inputfile1)
x.sort()
y = np.load(inputfile2)
y.sort()
plot_PDF(x, y,'Data from AF-cluster runs w/o templates',output_file2)

