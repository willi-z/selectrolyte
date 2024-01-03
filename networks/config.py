
import openpnm as op

class NetworkModelConfig:
    network: op.network.Network
    bounds: tuple[float, float, float]
    porosity: float
    def __init__(self, network: op.network.Network, bounds: tuple[float, float, float], porosity: float):
        self.network = network
        self.bounds = bounds
        self.porosity = porosity