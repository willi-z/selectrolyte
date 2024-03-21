import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from pathlib import Path
import json
import math

#"""
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'text.usetex': True,
    'pgf.rcfonts': False,
})
#"""

plot_config={
    "lines.color": "black",
    "patch.edgecolor": "black",
    "text.color": "black",
    
    "axes.linewidth": 1.5,
    "axes.facecolor": "white",

    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "figure.facecolor": "white",
    "figure.edgecolor": "black",

    # grid
    "axes.grid": True,
    "grid.color": "grey",
    "grid.linestyle": "dashed",
    
    # legend
    "legend.fancybox": False,
    "legend.edgecolor": "black",
    "legend.facecolor": "white",
    "legend.labelcolor": "black",
    "legend.framealpha": 0.8,
    
    "savefig.facecolor": "white",
    "savefig.edgecolor": "black",
    "savefig.transparent": True,
}

plt.rcParams.update(plot_config)
fig,ax = plt.subplots(figsize=(4,3))

num_pore = 3000

specimens = ["O1_50", "O2_40", "O4_40"]


with (Path.cwd() / f"data/data.json").open("r") as fp:
    data = json.load(fp)


inputDir = Path("/home/willi/Nextcloud/HTWK/share/selectrolyte/")
outputDir = inputDir
# process data
# Path.cwd() / "data/studies/err_over_ratio.json"
with (inputDir / "study_ratio.json").open("r") as fp:
    study = json.load(fp)

for specimen in specimens:
    Deff_exact = data[specimen]["conductivity"]
    xs = []
    ys = []
    min_distance = math.inf
    ratio_min_distance = -1

    ratios = study[specimen][str(num_pore)]
    for ratio, err in ratios.items():
        print(float(ratio))
        xs.append(float(ratio))
        err = np.array(err)
        rel_err = err / Deff_exact * 100 
        ys.append(np.mean(rel_err))
        # ys.append(rel_err[0])
        if abs(ys[-1]) < min_distance:
            min_distance = abs(ys[-1])
            ratio_min_distance = xs[-1]
    print(specimen, ":", ratio_min_distance)
    ax.plot(xs, ys, label=specimen)
# ax.set_xscale('log')
ax.set_xlabel(r"$r \; \left[ \; \right]$")
ax.set_ylabel(r"$rel. err \; \left[ \% \right]$")
#"""
ax.legend(loc='lower right', 
          fancybox=False, frameon=True)

fig.tight_layout(pad=0.1)

plt.savefig(outputDir / "study_ratio.png")