import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from pathlib import Path
import json

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

Deff_exact = 0.000005321864708362311
plt.rcParams.update(plot_config)
fig,ax = plt.subplots(figsize=(4,3))

num_pore = 1000

inputDir = Path("/home/willi/Nextcloud/HTWK/share/selectrolyte/bc_intersecting")
outputDir = inputDir
# process data
# Path.cwd() / "data/studies/err_over_ratio.json"
with (inputDir / "err_over_ratio.json").open("r") as fp:
    data = json.load(fp)

xs = []
ys = []

ratios = data[str(num_pore)]
for ratio, err in ratios.items():
    print(float(ratio))
    xs.append(float(ratio))
    err = np.array(err)
    rel_err = err / Deff_exact * 100 
    ys.append(np.mean(rel_err))
    # ys.append(rel_err[0])

print(len(xs), len(ys))

# ax.boxplot(x=ys, positions=xs, widths=0.01)
ax.plot(xs, ys)

ax.set_xlabel(r"$r \; \left[ \; \right]$")
ax.set_ylabel(r"$rel. err \; \left[ \% \right]$")
"""
ax.legend(loc='lower right', 
          fancybox=False, frameon=True)
"""

fig.tight_layout(pad=0)

#outputDir = Path().cwd() / "results"
plt.savefig(outputDir / (f"study_err_over_ratio_{num_pore}" + '.pdf'))
# plt.show()