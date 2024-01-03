from .ipore import IPoreGenerator, PoreNetworkConfig
from scipy.spatial import KDTree
import numpy as np
from .helpers import iso_2_coord


class RandomEquidistantSpacePoreGenerator(IPoreGenerator):
    def __init__(self)-> None:
        pass

    def generate(self, diameters: list[list[float]], bounds: tuple[float, float, float]) -> PoreNetworkConfig:
        drops = []
        coords = [iso_2_coord(np.random.rand(3), bounds)]
        diameters[::-1].sort()
        Dmax = diameters[0]

        for i in range(1, len(diameters)):
            tree = KDTree(
                coords,
            )
            Dthis = diameters[i]
            radius_min = (Dmax + Dthis)/2 *1.2
            # radius = box_length
            success = True
            for _ in range(100*i):
                success = True
                pos = iso_2_coord(np.random.rand(3), bounds)

                results = tree.query_ball_point(pos, radius_min)
                if len(results) == 0:
                    break
                for idx in results:
                    Dneigh = diameters[idx]
                    dist_squard = ((coords[idx] - pos)**2).sum()
                    if dist_squard < ((Dneigh + Dthis)/2 * 1.1)**2:
                        success = False
                        break
            if success:
                coords.append(pos)
            else:
                drops.append(i)

        Dpores = np.delete(diameters, drops)
        return PoreNetworkConfig(coords=coords, diameters=Dpores)