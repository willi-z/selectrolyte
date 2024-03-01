from pydantic import BaseModel


class PoreNetworkConfig(BaseModel):
    coords: list[list[float]]
    diameters: list[float]
    pore_ids: list[int]
    vol: float
    