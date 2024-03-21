import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


specimens = ["O1_50", "O2_40", "O4_40"]

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

fig,ax = plt.subplots(figsize=(4,3))

plt.rcParams.update(plot_config)

for specimen in specimens:
    xs = []
    ys = []
    distrib_file = './data/poresizes/dv_dlog/' + specimen + '.csv'
    with open(distrib_file, 'r') as csvfile:
        reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
        _header = reader.__next__()
        # print(header)
        for data in reader:
            xs.append(np.float64(data[0]))
            ys.append(np.float64(data[1]))

    ys = np.array(ys, dtype=np.float64)
    xs = np.array(xs, dtype=np.float64)
    ax.plot(xs, ys, label=specimen)

ax.set_xscale('log')
ax.set_xlabel(r"$w \; \left[ nm \right]$")
ax.set_ylabel(r"$dV/dlog(w) \; \left[ cm^3/g \right]$")
ax.legend(loc="upper right")
fig.tight_layout(pad=0.1)

plt.savefig('./data/graphics/gas_absorption.png')