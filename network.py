from networks.helpers import generate_network, update_network_connections  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import json
from pathlib import Path
import numpy as np


porosity = 0.54
num_pore = 1000
ratio = 0.001

distrib_file = './data/pore_distr_data.csv'
# np.random.seed(0)

conf = generate_network(
    CSVFileDiameterDistribution(distrib_file), # pore sizes in m
    num_pore,
    porosity,
    ratio,
    RandomPoreGenerator(),
    NearestConnsGenerator()
)

config_file = f"data/networks/rand_{porosity}.json"


with (Path.cwd() / config_file).open("w+") as fp:
    json.dump(conf.model_dump(), fp)

