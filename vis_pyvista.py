import pyvista
import matplotlib
import numpy as np
from pathlib import Path
import json
from networks.generators.config import NetworkConfig
from networks.helpers import model_from_network



# pc.plot(cmap='Reds')
boundary = pyvista.MultipleLines(points=[
    [0, 0, 0], 
    [0, 0, 1], 
    [0, 1, 1], 
    [0, 1, 0], 
    [0, 0, 0],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, 1],
    [1, 1, 0],
    [1, 0 ,0],
    [1, 0, 1],
    [0, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
    [1, 1, 0],
    [0, 1, 0]
])

pl = pyvista.Plotter()
pl.add_mesh(boundary)

with (Path.cwd() / "data/networks/random_0.54.json").open("r") as fp:
    content = json.load(fp)

conf = NetworkConfig(**content)
model = model_from_network(conf)
net = model.network
L = model.bounds[0]

DSpheres = net["pore.diameter"]
coords = net["pore.coords"]

norm = matplotlib.colors.Normalize(vmin=np.min(DSpheres),vmax=np.max(DSpheres))
cmap = matplotlib.cm.get_cmap("rainbow")

for i in range(len(DSpheres)):
    mesh = pyvista.Sphere(DSpheres[i]/2/L)
    mesh.translate((coords[i]/L), inplace=True)
    pl.add_mesh(mesh, color=cmap(norm(DSpheres[i])))

DThroats = net["throat.diameter"]
conns = net["throat.conns"]
for i in range(len(conns)):
    conn = conns[i]
    pstart = coords[conn[0]]/L
    pend = coords[conn[1]]/L
    mesh = pyvista.Tube(pstart, pend, radius=DThroats[i]/2/L)
    pl.add_mesh(mesh)
pl.show()