import numpy as np


def iso_2_coord(iso: np.ndarray, bounds: tuple[float, float, float]) -> np.ndarray:
    result = np.ones(3)
    result[0] = iso[0] * bounds[0]
    result[1] = iso[1] * bounds[1]
    result[2] = iso[2] * bounds[2]
    return result

def get_pore_vol(Dpores: np.ndarray):
    volumes = np.pi/6 * Dpores**3
    return volumes.sum()
