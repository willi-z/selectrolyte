from .config import ConnsNetworkConfig
from ..pores.config import PoreNetworkConfig

class IConnsGenerator:
    def generate(self, box_vol:float, porosity: float, pores: PoreNetworkConfig) -> ConnsNetworkConfig:
        raise NotImplementedError