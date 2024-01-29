import pyvista
import matplotlib
import numpy as np
from pathlib import Path
import json
from networks.generators.config import NetworkConfig
from networks.helpers import model_from_network
from networks.generators.pores.voxel import VoxelGrid



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

with (Path.cwd() / "data/networks/eq_0.54.json").open("r") as fp:
    content = json.load(fp)

conf = NetworkConfig(**content)
model = model_from_network(conf)
net = model.network
L = model.bounds[0]

DSpheres = net["pore.diameter"]
coords = net["pore.coords"]

norm = matplotlib.colors.Normalize(vmin=np.min(DSpheres),vmax=np.max(DSpheres))
cmap = matplotlib.cm.get_cmap("rainbow")



grid = VoxelGrid(conf.bounds, conf.pores.min_diameter)
drawn_indices = set()

for i in range(len(DSpheres)):
    removed_voxels = grid.remove_filled_voxels(coords[i], DSpheres[i]/2)
    drawn_indices.update(removed_voxels)

    mesh = pyvista.Sphere(DSpheres[i]/2/L)
    mesh.translate((coords[i]/L), inplace=True)
    pl.add_mesh(mesh, opacity= 0.5, color=cmap(norm(DSpheres[i])))

voxel_length = (
    grid.voxel_sizes[0]/L,
    grid.voxel_sizes[1]/L,
    grid.voxel_sizes[2]/L,
)
for index in drawn_indices:
    center = grid.get_voxel_center(index)
    center = (
        center[0] / L,
        center[1] / L,
        center[2] / L
    )
    mesh = pyvista.Cube(center, voxel_length[0], voxel_length[1], voxel_length[2])
    pl.add_mesh(mesh, opacity= 0.8)


pl.show()