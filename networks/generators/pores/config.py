from pydantic import BaseModel


class PoreNetworkConfig(BaseModel):
    coords: list[list[float]]
    diameters: list[float]
    extra_volume: float = 0.0
    min_diameter: float = 0.0
    