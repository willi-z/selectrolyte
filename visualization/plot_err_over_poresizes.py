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
    #"axes.grid": True,
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
fig,ax = plt.subplots(figsize=(4,3), dpi=600)

num_pore = 3000

specimens = ["O1_50", "O2_40", "O4_40"]


with (Path.cwd() / f"data/data.json").open("r") as fp:
    data = json.load(fp)

ratio_throat_pores = {
    specimens[0]: 0.1,
    specimens[1]: 0.007,
    specimens[2]: 0.52
}

# process data
inputDir = Path("/home/willi/Nextcloud/HTWK/share/selectrolyte")
with (inputDir / "study_poresize.json").open("r") as fp:
    study = json.load(fp)

#for specimen in specimens:
specimen_id = 2
specimen = specimens[specimen_id]
if True:
    Deff_exact = data[specimen]["conductivity"]
    ratio = ratio_throat_pores[specimen]
    xs = []
    ys = []
    values = study[specimen]
    for num_pore, ratios in values.items():
        if len(ratios) == 0:
            continue
        #if int(num_pore) < 1000:
        #    continue    
        xs.append(int(num_pore))
        print(num_pore)
        err = np.array(ratios[str(ratio)])
        rel_err = err / Deff_exact * 100 
        ys.append(rel_err)

    ax.boxplot(x=ys, positions=xs, widths=100, medianprops=dict(color=list(plt.rcParams['axes.prop_cycle'])[specimen_id]["color"], alpha=0.7),)

ax.set_xlabel(r"$n_p \; \left[ \; \right]$")
ax.set_ylabel(r"$rel. err \; \left[ \% \right]$")
"""
ax.legend(loc='lower right', 
          fancybox=False, frameon=True)
"""

fig.tight_layout(pad=0.1)

outputDir = Path().cwd() / "results"
plt.savefig(inputDir / (f"study_poresizes_{specimen}.png"))
# plt.show()