from pathlib import Path
import json
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

specimen = "O1_50_3000_0.1"

with (Path.cwd() / "data/networks" / (specimen + ".json")).open("r") as fp:
    content = json.load(fp)

conf = NetworkConfig(**content)
model = model_from_network(conf, BC_Scale=1.0)
net = model.network

plt.clf()
plt.cla()
fig, ax = plt.subplots(figsize=(4,3))
plt.rcParams.update(plot_config)
ax.hist(net['pore.diameter']* 1e9, bins=25, edgecolor='k')
ax.set_xlabel(r"$D_p \; \left[ nm \right]$")
ax.set_ylabel(r"$n_p \; \left[ \; \right]$")
fig.tight_layout(pad=0.1)
plt.savefig("./data/graphics/" + specimen + '_pore.png')
plt.clf()
plt.cla()


fig, ax = plt.subplots(figsize=(4,3))
plt.rcParams.update(plot_config)
ax.hist(net['throat.diameter'] * 1e9, bins=25, edgecolor='k')
ax.set_xlabel(r"$D_t \; \left[ nm \right]$")
ax.set_ylabel(r"$n_t \; \left[ \; \right]$")
fig.tight_layout(pad=0.1)
plt.savefig("./data/graphics/" + specimen + '_throat.png')