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

ratio = 0.06

# process data
with (Path.cwd() / "data/studies/err_over_ratio.json").open("r") as fp:
    data = json.load(fp)

xs = []
ys = []
for num_pore, ratios in data.items():
    if len(ratios) == 0:
        continue
    xs.append(int(num_pore))
    print(num_pore)
    err = np.array(ratios[str(ratio)])
    rel_err = err / Deff_exact * 100 
    ys.append(rel_err)


ax.boxplot(x=ys, positions=xs, widths=20)

ax.set_xlabel(r"$n_p \; \left[ \; \right]$")
ax.set_ylabel(r"$rel. err \; \left[ \% \right]$")
"""
ax.legend(loc='lower right', 
          fancybox=False, frameon=True)
"""

fig.tight_layout(pad=0)

outputDir = Path().cwd() / "results"
plt.savefig(outputDir / (f"study_err_over_poresizes_{ratio}" + '.pdf'))
# plt.show()