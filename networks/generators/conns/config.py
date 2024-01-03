from pydantic import BaseModel


class ConnsNetworkConfig(BaseModel):
    conns: list[tuple[int, int]]
    diameters: list[float]