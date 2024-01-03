from .config import PoreNetworkConfig

class IPoreGenerator:
    def generate(self, diameters: list[list[float]], bounds: tuple[float, float, float]) -> PoreNetworkConfig:
        raise NotImplementedError