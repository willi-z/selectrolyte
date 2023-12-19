import csv
import numpy as np
from scipy.spatial import KDTree
import openpnm as op

num_pores = 200
porosity = 0.54
eps = 0.05

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

def iso_2_coord(iso):
    return iso*box_length


coords = [iso_2_coord(np.random.rand(3))]

drops = []
print(box_length)
print(Dmax)

for i in range(1, len(Dpores)):
    tree = KDTree(
        coords,
    )
    radius_min = (Dmax + Dpores[i])/2 *1.2
    # radius = box_length
    success = True
    for _ in range(100*i):
        success = True
        pos = iso_2_coord(np.random.rand(3))

        results = tree.query_ball_point(pos, radius_min)
        if len(results) == 0:
            break
        for idx in results:
            Dneigh = Dpores[idx]
            dist_squard = ((coords[idx] - pos)**2).sum()
            if dist_squard < ((Dneigh + Dpores[i])/2 * 1.1)**2:
                success = False
                break
    if success:
        coords.append(pos)
    else:
        drops.append(i)

Dpores = np.delete(Dpores, drops)

print(len(coords))
print(len(Dpores))

tree = KDTree(
        coords
)
tubes = {}
for i in range(len(coords)):
    pos = coords[i]
    
    _, results = tree.query(pos, 3)
    neighbours = set(results)
    neighbours.remove(i)
    for idx in neighbours:
        Dneigh = Dpores[idx]
        DTube = min(Dneigh, Dpores[i])
        tubes[frozenset([idx, i])]=DTube

    radius_max = (Dpores[i] + Dmax)/2 *1.6
    results = tree.query_ball_point(pos, radius_max)
    neighbours = set(results)
    neighbours.remove(i)
    print(len(neighbours))
    if neighbours is None:
        continue
    for idx in neighbours:
        Dneigh = Dpores[idx]
        distance = np.sqrt(((coords[idx] - pos)**2).sum())
        distance_max = (Dneigh + Dpores[i])/2 * 1.5
        if distance > distance_max:
            continue
        DTube = min(Dneigh, Dpores[i])
        tubes[frozenset([idx, i])]=DTube


connections = []
Dthroats = []
for idxs, dia in tubes.items():
    conns = list(idxs)
    connections.append(conns)
    distance = np.sqrt(((coords[conns[0]] - coords[conns[1]])**2).sum())
    Dthroats.append(dia*(1-dia))


print(len(connections))

pore_vol = get_pore_vol()
print(pore_vol/box_vol)

net = op.network.Network(coords=coords, conns=connections)
# print(len(coords) == len(Dpores), len(coords), len(Dpores))
net['pore.diameter'] = Dpores
net['throat.diameter'] = Dthroats

net.add_model(propname='pore.volume',
             model=op.models.geometry.pore_volume.sphere)
net.add_model(propname='throat.length',
             model=op.models.geometry.throat_length.spheres_and_cylinders)
net.add_model(propname='throat.total_volume',
             model=op.models.geometry.throat_volume.cylinder)
net.add_model(propname='throat.lens_volume', 
            model=op.models.geometry.throat_volume.lens)
net.add_model(propname='throat.volume', 
             model=op.models.misc.difference,
             props=['throat.total_volume', 'throat.lens_volume'])


net.add_model(
    propname='throat.hydraulic_size_factors', 
    model=op.models.geometry.hydraulic_size_factors.spheres_and_cylinders
)

net.add_model(
    propname='throat.diffusive_size_factors', 
    model=op.models.geometry.diffusive_size_factors.spheres_and_cylinders
)
net.regenerate_models()
xPore = np.array(coords)[:, 0]
net["pore.left"] = np.isclose(xPore, 0.0, atol=box_length*eps)
net["pore.right"] = np.isclose(xPore, box_length, atol=box_length*eps)

print("left: ", np.where(net['pore.left'])[0])
print("right: ", np.where(net['pore.right'])[0])

Vol_void = np.sum(net['pore.volume'])+np.sum(net['throat.volume'])
Vol_bulk = box_vol
porosity = Vol_void / Vol_bulk
print(f'The value of Porosity is: {porosity * 100:.2f}%')

A = box_length * box_length
L = box_length