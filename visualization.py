import plotly.graph_objects as go
from networks.randnet import net
import numpy as np

fig = go.Figure()
coords = net["pore.coords"]
xs = np.zeros(len(coords))
ys = np.zeros(len(coords))
zs = np.zeros(len(coords))

sizes = np.array(net["pore.diameter"])
max_diameter = max(np.array(sizes))
norm_sizes = sizes / max(np.array(sizes))

for i in range(len(coords)):
    coord = coords[i]
    xs[i] = coord[0]
    ys[i] = coord[1]
    zs[i] = coord[2]

#"""
fig.add_trace(go.Scatter3d(
    x = xs,
    y = ys,
    z = zs,
    mode='markers',
    marker=dict(
        color=sizes,
        colorscale='bluered',
        opacity=1,
        size=norm_sizes * 8.0,
    
        colorbar=dict(
            title="Pore Diameter [m]",
            titleside="top",
            tickmode="array",
            #tickvals=[2, 25, 50, 75, 100],
            #labelalias={100: "Hot", 50: "Mild", 2: "Cold"},
            #ticks="outside"
        )
    ),
    showlegend=False,
))
#"""

"""
connections = np.array(net["throat.conns"])
sizes = np.array(net["throat.diameter"])
norm_sizes = sizes / max_diameter
for i in range(len(connections)): # len(connections)
    con = connections[i]
    # print(con)
    fig.add_trace(go.Scatter3d(
    x = [xs[con[0]], xs[con[1]]],
    y = [ys[con[0]], ys[con[1]]],
    z = [zs[con[0]], zs[con[1]]],
    mode='lines',
    line=dict(
        colorscale="bluered",
        width= norm_sizes[i] * 8.0
    ),
    showlegend=False,
))
#"""

# fig.show()
fig.write_html("results/network_pores_only.html")

