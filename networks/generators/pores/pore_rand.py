from .ipore import IPoreGenerator, PoreNetworkConfig
from scipy.spatial import KDTree
import numpy as np
from .helpers import iso_2_coord


class RandomPoreGenerator(IPoreGenerator):
    def __init__(self) -> None:
        pass

    def generate(
        self, diameters: list[list[float]], bounds: tuple[float, float, float]
    ) -> PoreNetworkConfig:
        diameters[::-1].sort()
        Dmax = diameters[0]
        Dmin = diameters[-1]

        drops = []
        # np.random.seed(0)
        coords = [iso_2_coord(np.random.rand(3), ((0, 0, 0), bounds))]

        # grid = VoxelGrid(bounds, Dmin)
        # grid.remove_filled_voxels(coords[0], diameters[0])

        for i in range(1, len(diameters)):
            tree = KDTree(
                coords,
            )
            Dthis = diameters[i]
            radius_min = (Dmax + Dthis) / 2 * 1.2
            # radius = box_length
            success = False
            pos = None
            counter = 0
            # print(i, ": ", end="")
            while counter < 1000 and not success:  # and counter < 1000:
                # print(len(possible_voxel_indices), end="\r")
                counter += 1
                success = True
                pos = iso_2_coord(np.random.rand(3), ((0, 0, 0), bounds))

                results = tree.query_ball_point(pos, radius_min)
                for idx in results:
                    Dneigh = diameters[idx]
                    dist_squard = ((coords[idx] - pos) ** 2).sum()
                    if dist_squard < ((Dneigh + Dthis) / 2 * 1.1) ** 2:
                        success = False
                        break
            # print(counter, "(success:", success, ")")
            if success and pos is not None:
                # assert pos is not None
                # grid.remove_filled_voxels(pos, diameters[i]/2)
                coords.append(pos)
            else:
                drops.append(i)

        # print("droped dias (larger index = smaller):", drops)
        Dmin = diameters[-1]
        Dpores = np.delete(diameters, drops)[: len(coords)]
        return PoreNetworkConfig(coords=coords, diameters=Dpores, min_diameter=Dmin)
