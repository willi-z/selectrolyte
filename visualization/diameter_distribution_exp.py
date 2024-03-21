from pathlib import Path
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


import imp
networks = imp.load_source('networks', str(Path.cwd() / "networks/__init__.py"))

from networks.generators.config import NetworkConfig
from networks.helpers import model_from_network


plot_config={
    "figure.dpi": 300,
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
    "axes.grid": False,
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

specimens = ["O1_50", "O2_40", "O4_40"]
specimen = specimens[0]
distrib_file = './data/poresizes/comulative_pore_volume/' + specimen + '.csv'
pore_no = 3000

xs = []
ys = []

with open(distrib_file, 'r') as csvfile:
    reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
    _header = reader.__next__()
    # print(header)
    for data in reader:
        xs.append(np.float64(data[0]))
        ys.append(np.float64(data[1]))
        # ys.append(np.float64(data[0]))

    ys = np.array(ys, dtype=np.float64) - ys[0]
    ys = ys / ys[-1] # normalized
    xs = np.array(xs, dtype=np.float64) * 1e-1 # * 1e-10 # A to m


seeds = np.random.rand(pore_no)
Dpores = np.interp(seeds, ys, xs)


fig, ax = plt.subplots(figsize=(4,3))
plt.rcParams.update(plot_config)
ax.hist(Dpores, bins=25, edgecolor='k')
ax.set_xlabel(r"$D_p \; \left[ nm \right]$")
ax.set_ylabel(r"$n_p \; \left[ \; \right]$")
fig.tight_layout(pad=0.1)
plt.savefig("./data/graphics/" + f"{specimen}_{pore_no}" + '_pore_exp.png')
