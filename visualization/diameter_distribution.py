from pathlib import Path
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


import imp
networks = imp.load_source('networks', str(Path.cwd() / "networks/__init__.py"))

from networks.generators.config import NetworkConfig
from networks.helpers import model_from_network


matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'text.usetex': True,
    'pgf.rcfonts': False,
})

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


with (Path.cwd() / "data/networks/rand_0.54.json").open("r") as fp:
    content = json.load(fp)

conf = NetworkConfig(**content)
model = model_from_network(conf, BC_Scale=1.0)
net = model.network

fig, ax = plt.subplots(figsize=(4,3))
ax.hist(net['pore.diameter']* 1e10, bins=25, edgecolor='k')
ax.set_xlabel(r"$D_p \; \left[ \buildrel _\circ \over {\mathrm{A}} \right]$")
ax.set_ylabel(r"$n_p \; \left[ \; \right]$")
fig.tight_layout(pad=0)
plt.savefig(f"./data/diameter_distribution_pore" + '.pdf')
plt.clf()

fig, ax = plt.subplots(figsize=(4,3))
ax.hist(net['throat.diameter'] * 1e10, bins=25, edgecolor='k')
ax.set_xlabel(r"$D_t \; \left[ \buildrel _\circ \over {\mathrm{A}} \right]$")
ax.set_ylabel(r"$n_t \; \left[ \; \right]$")
fig.tight_layout(pad=0)
plt.savefig(f"./data/diameter_distribution_throat" + '.pdf')