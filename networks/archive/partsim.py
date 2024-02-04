import open3d as o3d
import numpy as np
import csv
import numpy as np
from numpy import linalg as LA
from scipy.spatial import KDTree
import openpnm as op
import matplotlib

"""
https://www.youtube.com/watch?v=rSKMYc1CQHE
https://stackoverflow.com/questions/65774814/adding-new-points-to-point-cloud-in-real-time-open3d
"""

# Global settings.
dt = 1e-6 # to add new points each dt secs.
t_total = 10 # total time to run this script.
n_new = 10 # number of points that will be added each iteration.
num_pores = 5
porosity = 0.4
eps = 0.05
colli_damp = 0.2
density_target = 1e-5
p_multiplier = 1e-4

xs = [] # pore size
ys = [] # comulative volume
with open('./data/pore_distr_data.csv', 'r') as csvfile:
    reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
    header = reader.__next__()
    print(header)
    for data in reader:
        xs.append(np.float64(data[0]))
        ys.append(np.float64(data[1]))

ys = np.array(ys, dtype=np.float64) - ys[0]
ys = ys / ys[-1] # normalized
xs = np.array(xs, dtype=np.float64) * 1e-10 # A to m

np.random.seed(0)

seeds = np.random.rand(num_pores)

Dpores = np.interp(seeds, ys, xs)
# sort Dpores form largest to smallest
Dpores[::-1].sort()

Dmax = Dpores[0]

def get_pore_vol():
    volumes = np.pi/6 * Dpores**3
    return volumes.sum()

pore_vol = get_pore_vol()

box_vol = pore_vol / porosity
box_length = np.power(box_vol,1/3.)
print(pore_vol / (box_length**3))

def iso_2_coord(iso):
    return iso*box_length

def coord_2_iso(coord):
    return coord/box_length


coords = []
coords_predicted =  []
for i in range(len(Dpores)):
    coords.append(iso_2_coord(np.random.rand(3)))
    coords_predicted.append(coords[i].copy())

velocities = np.zeros((num_pores, 3), np.float32)
densities = np.zeros((num_pores,2), np.float32)

def resolve_collisons(idx):
    coord  = coords[idx]
    for i in range(3):
        if coord[i] < 0.0:
            coord[i] = 0.0
            velocities[idx][i] *= -1 * colli_damp
        elif coord[i] > box_length:
            coord[i] = box_length
            velocities[idx][i] *= -1 * colli_damp

def spiky_kernel_pow2(dist: float, radius: float)->float:
    if dist > radius:
        return 0
    scale = 15 / (2 * np.pi * np.power(radius, 5))
    v = radius - dist
    return v * v * scale

def spiky_kernel_pow3(dist: float, radius: float)->float:
    if dist > radius:
        return 0
    scale = 15 / (2 * np.pi * np.power(radius, 6))
    v = radius - dist
    return v * v * v * scale

def spiky_kernel_pow2_deriv(dist: float, radius: float)->float:
    if dist > radius:
        return 0
    scale = 15 / (np.pi * np.power(radius, 5))
    v = radius - dist
    return -v * scale

def spiky_kernel_pow3_deriv(dist: float, radius: float)->float:
    if dist > radius:
        return 0
    scale = 45 / (np.pi * np.power(radius, 6))
    v = radius - dist
    return -v * v * scale


def density_kernel(dist: float, radius: float):
    return spiky_kernel_pow2(dist, radius)

def near_density_kernel(dist: float, radius: float):
    return spiky_kernel_pow3(dist, radius)

def density_kernel_deriv(dist: float, radius: float):
    return spiky_kernel_pow2_deriv(dist, radius)

def near_density_kernel_deriv(dist: float, radius: float):
    return spiky_kernel_pow3_deriv(dist, radius)

def calc_mass(idx):
    volume = np.pi / 6 * Dpores[idx]
    return volume * density_target

def calc_density(tree, idx):
    global densities
    density = 0.0
    near_density = 0.0
    smooth_radius = box_length # Dmax # Dpores[idx]
    pos = coords_predicted[idx]
    neighbours = tree.query_ball_point(pos, smooth_radius)
    for nidx in neighbours:
        npos = coords_predicted[nidx]
        mass = calc_mass(nidx)
        offset = npos - pos
        dist = LA.norm(offset)
        density += density_kernel(dist, smooth_radius) * mass
        near_density += near_density_kernel(dist, smooth_radius) * mass
    densities[idx][0], densities[idx][1] = density, near_density

def pressure_from_density(density: float):
    density_err = density - density_target
    pressure = density_err * p_multiplier
    return pressure

def near_pressure_from_density(near_density: float):
    pressure = near_density * p_multiplier
    return pressure

def calc_pressure_force(tree, idx):
    global densities, velocities
    pForce = np.zeros(3)
    density, near_density = tuple(densities[idx])
    pressure, near_pressure = pressure_from_density(density), near_pressure_from_density(near_density)
    smooth_radius = box_length# Dmax # Dpores[idx]
    pos = coords_predicted[idx]
    neighbours = tree.query_ball_point(pos, smooth_radius)
    for nidx in neighbours:
        if nidx == idx:
            continue

        offset = coords_predicted[nidx] - pos
        dist = LA.norm(offset)

        direction = np.random.rand(3)
        if dist == 0.0:
            direction = direction / LA.norm(direction)
        else:
            direction = offset / dist
            
        density_neighbour, near_density_neighbour = densities[nidx]
        pressure_neighbour, near_pressure_neighbour = pressure_from_density(density_neighbour), near_pressure_from_density(near_density_neighbour)

        pressure_share = (pressure + pressure_neighbour) / 2
        near_pressure_share = (near_pressure + near_pressure_neighbour) / 2

        pf1 = direction * density_kernel_deriv(dist, smooth_radius) * pressure_share / density_neighbour
        pf2 = direction * near_density_kernel_deriv(dist, smooth_radius) * near_pressure_share / near_density_neighbour

        pForce += pf1
        pForce += pf2
    accelaration = pForce / density
    velocities[idx] += accelaration * dt

def update_coords(coords):
    tree = KDTree(coords_predicted)
    for i in range(num_pores):
        # velocities[i] += np.array([0.0,-1.0,0.0]) * (0.1*box_length) * dt
        calc_density(tree, i)

    for i in range(num_pores):
        calc_pressure_force(tree, i)

    for i in range(num_pores):
        coords[i] += velocities[i] *dt
        resolve_collisons(i)
        coords_predicted[i] = coords[i] + velocities[i] / 120.0
    return coords




# visualization
norm = matplotlib.colors.Normalize(vmin=np.min(Dpores),vmax=np.max(Dpores))
cmap = matplotlib.colormaps["rainbow"]

Mpores = []

for i in range(len(Dpores)):
    mesh = o3d.geometry.TriangleMesh.create_sphere(radius=Dpores[i]/2/box_length)
    color = np.array(cmap(norm(Dpores[i])))
    rgb = color[:3]
    mesh.paint_uniform_color(rgb)
    Mpores.append(mesh)

#---
# 1st, using extend. Run non-blocking visualization.

# create visualizer and window.
vis = o3d.visualization.Visualizer()
vis.create_window(window_name="Pore Network", height=600, width=600)

def generate_border():
    points = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ]
    lines = [
        [0, 1],
        [0, 2],
        [1, 3],
        [2, 3],
        [4, 5],
        [4, 6],
        [5, 7],
        [6, 7],
        [0, 4],
        [1, 5],
        [2, 6],
        [3, 7],
    ]
    colors = [[0.2, 0.2, 0.2] for i in range(len(lines))]
    border = o3d.geometry.LineSet(
        points=o3d.utility.Vector3dVector(points),
        lines=o3d.utility.Vector2iVector(lines),
    )
    border.colors = o3d.utility.Vector3dVector(colors)
    return border
vis.add_geometry(generate_border())

# include it in the visualizer before non-blocking visualization.
# vis.add_geometry(pcd)

for mesh in Mpores:
    vis.add_geometry(mesh)

MCpores = np.zeros((num_pores, 3), np.float32)
def update_mesh(coords):
    for i in range(num_pores):
        mesh = Mpores[i]
        coord_new = coords[i]
        coord_old = MCpores[i]
        dcoord = (coord_new - coord_old)/box_length
        MCpores[i] = coord_new
        mesh.translate(dcoord)
        vis.update_geometry(mesh)
update_mesh(coords)
    


for _ in range(50000):
    coords = update_coords(coords)
    update_mesh(coords)
    vis.poll_events()
    vis.update_renderer()
#vis.poll_events()
# vis.update_renderer()
vis.run()