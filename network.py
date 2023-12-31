from networks.helpers import generate_network, update_network_connections  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import json
from pathlib import Path


num_pores = 100
porosity = 0.54
ratio_throat_pores = 0.02

distrib_file = './data/pore_distr_data.csv'
config_file = f"data/networks/rand_{porosity}.json"

# """"
conf = generate_network(
    CSVFileDiameterDistribution(distrib_file),
    num_pores,
    porosity,
    ratio_throat_pores,
    RandomPoreGenerator(),
    NearestConnsGenerator()
)


"""

with (Path.cwd() / config_file).open("r") as fp:
    content = json.load(fp)

conf = update_network_connections(
    NetworkConfig(**content),
    NearestConnsGenerator()
)
# """

print("num. pores:", len(conf.pores.diameters))
print("num. conncetions: ",len(conf.conns.conns))
with (Path.cwd() / config_file).open("w+") as fp:
    json.dump(conf.model_dump(), fp)

"""
A = conf.bounds[1] * conf.bounds[2]
L = conf.bounds[0]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[12, 4])
ax1.hist(net['pore.diameter'], bins=40, edgecolor='k')
ax1.set_title('Pore Diameter')
ax2.hist(net['throat.diameter'], bins=25, edgecolor='k')
ax2.set_title('Throat Diameter')
plt.savefig("results/network_distr.png")
"""