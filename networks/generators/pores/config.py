from pydantic import BaseModel


class PoreNetworkConfig(BaseModel):
    coords: list[list[float]]
    diameters: list[float]
    min_diameter: float = 0.0
    