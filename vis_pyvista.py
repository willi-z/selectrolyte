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

with (Path.cwd() / "data/networks/rand_0.54.json").open("r") as fp:
    content = json.load(fp)

conf = NetworkConfig(**content)
model = model_from_network(conf, BC_Scale=1.0)
net = model.network
L = conf.bounds[0]


DSpheres = conf.pores.diameters
coords = np.array(conf.pores.coords)

norm = matplotlib.colors.Normalize(vmin=np.min(DSpheres),vmax=np.max(DSpheres))
cmap = matplotlib.colormaps["rainbow"]



for i in range(len(DSpheres)): # len(DSpheres)
    opacity = 0.2
    if net["pore.left"][i] or net["pore.right"][i]:
        opacity = 1.0
    radius = DSpheres[i]/2/L
    mesh = pyvista.Sphere(radius)
    mesh.translate((coords[i]/L), inplace=True)
    pl.add_mesh(mesh, opacity=opacity, color=cmap(norm(DSpheres[i])))


# if conf.conns is None:
DThroats = conf.conns.diameters # ["throat.diameter"]
conns = conf.conns.conns # ["throat.conns"]
for i in range(len(conns)):
    conn = conns[i]
    pstart = coords[conn[0]]/L
    pend = coords[conn[1]]/L
    mesh = pyvista.Tube(pstart, pend, radius=DThroats[i]/2/L)
    pl.add_mesh(mesh)

# print(L)
# print(max(DThroats))

print("Finished!")
pl.show()