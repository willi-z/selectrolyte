from networks.helpers import generate_network, update_network_connections  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import json
from pathlib import Path


porosity = 0.54

config_file = f"data/networks/rand_{porosity}.json"


with (Path.cwd() / config_file).open("r") as fp:
    content = json.load(fp)

