from .config import ConnsNetworkConfig
from ..pores.config import PoreNetworkConfig

class IConnsGenerator:
    def generate(self, bounds:tuple[float, float, float], porosity: float, pores: PoreNetworkConfig) -> ConnsNetworkConfig:
        raise NotImplementedError