import pyvista
import matplotlib
import numpy as np
from pathlib import Path
import json
import imp
networks = imp.load_source('networks', str(Path.cwd() / "networks/__init__.py"))
from networks.generators.config import NetworkConfig
from networks.helpers import model_from_network



specimens = ["O1_50", "O2_40", "O4_40"]
specimen = specimens[0]
porosity_str = '0.50'
num_pore_str = '1000'
ratio_str = '0.02'


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

pl = pyvista.Plotter(window_size=(1200,1200))
# pl.view_isometric()
# pl.camera.zoom(0.25)
pl.add_mesh(boundary, color="k",line_width=2)
show_only_boundary = True

with (Path.cwd() / f"data/networks/{specimen}_{num_pore_str}_{porosity_str}_{ratio_str}.json").open("r") as fp:
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
    opacity = 1.0
    if show_only_boundary:
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
    opacity = 1.0
    if show_only_boundary:
        opacity = 0.2
        if net["pore.left"][conn[0]] or net["pore.right"][conn[0]]:
            opacity = opacity + 0.4
        if net["pore.left"][conn[1]] or net["pore.right"][conn[1]]:
            opacity = opacity + 0.4
    pl.add_mesh(mesh, opacity=opacity)

# print(L)
# print(max(DThroats))

print("Finished!")
pl.show()