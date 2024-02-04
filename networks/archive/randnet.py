import csv
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt


xs = [] # pore size
ys = [] # comulative volume
with open('../data/pore_distr_data.csv', 'r') as csvfile:
    reader = csv.reader(csvfile,quoting=csv.QUOTE_MINIMAL)
    header = reader.__next__()
    print(header)
    for data in reader:
        xs.append(np.float64(data[0]))
        ys.append(np.float64(data[1]))

ys = np.array(ys, dtype=np.float64) - ys[0]
ys = ys / ys[-1] # normalized
xs = np.array(xs, dtype=np.float64) # * 1e-10 # A to m

class GDistribution:
    def __init__(self) -> None:
        pass

    def ppf(self, seeds):
        return np.interp(seeds, ys, xs)

np.random.seed(0)
#shape = [1,1,1]
#spacing = 1.0
#net = op.network.DelaunayVoronoiDual(shape=shape, points=1000)

spacing, shape = 0.143e-7 * 1e10, [5,5,3]
# spacing, shape = 0.143e-7, [10,10,10]
# shape = [50,50,50]
# spacing, shape = 0.094e-6, [30,30,30]

net = op.network.Cubic(shape=shape, spacing=spacing)

f_random = op.models.geometry.pore_seed.random
net.add_model(propname='pore.seed',
             model=f_random,
             num_range=[0.0, 1.0])

dst = GDistribution()
net.add_model(propname='pore.diameter',
             model=op.models.geometry.pore_size.generic_distribution,
             func=dst,
             seeds='pore.seed')

net.add_model(propname='throat.diameter_1', 
             model=op.models.misc.from_neighbor_pores,
             prop='pore.diameter',
             mode='min')
net.add_model(propname='throat.seed', 
             model=op.models.misc.from_neighbor_pores,
             prop='pore.seed',
             mode='min')
net.add_model(propname='throat.diameter',
             model=op.models.misc.scaled,
             prop='throat.diameter_1',  # This could also be 'throat.diameter_2'
             factor=0.6,  # This could be 1.0 if no scaling is desired
)


diameters = net["pore.diameter"]
connections = net["throat.conns"]
drop = []
for i in range(len(connections)):
    connection = connections[i]
    d0 = diameters[connection[0]]
    d1 = diameters[connection[1]]

    if (d0 + d1)/2 < spacing *1.1:
        drop.append(i)

op.topotools.trim(network=net, throats=drop)
print

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

A = (shape[1] * shape[2])*(spacing**2)
L = shape[0]*spacing


Vol_void = np.sum(net['pore.volume'])+np.sum(net['throat.volume'])
Vol_bulk = shape[0] * spacing * shape[1] * spacing * shape[2] * spacing
porosity = Vol_void / Vol_bulk
print(f'The value of Porosity is: {porosity * 100:.2f}%')
# """
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[12, 4])
ax1.hist(net['pore.diameter'], bins=40, edgecolor='k')
ax1.set_title('Pore Diameter')
ax2.hist(net['throat.diameter'], bins=25, edgecolor='k')
ax2.set_title('Throat Diameter')
plt.savefig("results/network_distr.png")
# plt.show()
#"""