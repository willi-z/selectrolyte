from networks.helpers import generate_network, update_network_connections  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import json
from pathlib import Path
import numpy as np



specimens = ["O1_50", "O2_40", "O4_40"]
specimen = specimens[0]
porosity_str = '0.50'
num_pore_str = '1000'
ratio_str = '0.02'

porosity = float(porosity_str)
num_pore = int(num_pore_str)
ratio = float(ratio_str)
print(ratio)

distrib_file = './data/poresizes/comulative_pore_volume/' + specimen + '.csv'
# np.random.seed(0)

conf = generate_network(
    CSVFileDiameterDistribution(distrib_file), # pore sizes in m
    num_pore,
    porosity,
    ratio,
    RandomPoreGenerator(),
    NearestConnsGenerator()
)

config_file = f"data/networks/{specimen}_{num_pore_str}_{porosity_str}_{ratio_str}.json"


with (Path.cwd() / config_file).open("w+") as fp:
    json.dump(conf.model_dump(), fp)

