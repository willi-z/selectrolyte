import numpy as np


def iso_2_coord(
    iso: np.ndarray, 
    bounds: tuple[
        tuple[float, float, float],
        tuple[float, float, float]
    ]
    ) -> np.ndarray:
    bounds_min, bounds_max = bounds
    bounds_dx = (
        bounds_max[0] - bounds_min[0],
        bounds_max[1] - bounds_min[1],
        bounds_max[2] - bounds_min[2],
    )
    result = np.zeros(3)
    result[0] = iso[0] * bounds_dx[0] + bounds_min[0]
    result[1] = iso[1] * bounds_dx[1] + bounds_min[1]
    result[2] = iso[2] * bounds_dx[2] + bounds_min[2]
    return result

def get_pore_vol(Dpores: np.ndarray):
    volumes = np.pi/6 * Dpores**3
    return volumes.sum()
