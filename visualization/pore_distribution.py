import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


specimens = ["O1_50", "O2_40", "O4_40"]
specimen = specimens[0]
distrib_file = './data/poresizes/comulative_pore_volume/' + specimen + '.csv'

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


xs = []
ys = []

with open(distrib_file, 'r') as csvfile:
    reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
    _header = reader.__next__()
    # print(header)
    for data in reader:
        xs.append(np.float64(data[0]))
        ys.append(np.float64(data[1]))

ys = np.array(ys, dtype=np.float64) - ys[0]
ys = ys / ys[-1] # normalized
xs = np.array(xs, dtype=np.float64) # * 1e-10 # A to m


fig,ax = plt.subplots(figsize=(4,3))
ax.plot(ys, xs)
ax.set_ylabel(r"$D_p \; \left[ \buildrel _\circ \over {\mathrm{A}} \right]$")
ax.set_xlabel(r"$p \; \left[ \; \right]$")
fig.tight_layout(pad=0)

plt.savefig(f"./data/pore_distribution_" + specimen + '.pdf')