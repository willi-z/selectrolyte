from config import NetworkConfig
from ..distribution.idistribution import IDiameterDistribution

class NetworkGenerator:
    def __init__(self, diameter_generator: IDiameterDistribution):
        pass

    def generate(porosity: float, pore_generator) -> NetworkConfig:
        raise NotImplementedError
