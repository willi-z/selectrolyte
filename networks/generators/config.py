from pydantic import BaseModel
from .pores.config import PoreNetworkConfig
from .conns.config import ConnsNetworkConfig


class NetworkConfig(BaseModel):
    pores: PoreNetworkConfig
    conns: ConnsNetworkConfig
    bounds: tuple[float, float, float]
    porosity: float