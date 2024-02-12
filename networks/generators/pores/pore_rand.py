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
        coords = [] # iso_2_coord(np.random.rand(3), ((0, 0, 0), bounds))
 
        # grid = VoxelGrid(bounds, Dmin)
        # grid.remove_filled_voxels(coords[0], diameters[0])
        cidx = 0

        while cidx < len(diameters):
            
            Dthis = diameters[cidx]
            radius_min = (Dmax + Dthis) / 2 * 1.2
            # radius = box_length
            success = False
            pos = None
            shadow_pores = [] # in case an overhang occures

            counter = 0

            tree = None
            if len(coords) != 0:
                tree = KDTree(
                    coords,
                )
            while counter < 1000 and not success:
                counter += 1
                success = True
                iso_coord = np.random.rand(3)
                pos = iso_2_coord(iso_coord, ((0, 0, 0), bounds))

                if tree is not None:
                    results = tree.query_ball_point(pos, radius_min)
                    for idx in results:
                        Dneigh = diameters[idx]
                        dist_squard = ((coords[idx] - pos) ** 2).sum()
                        if dist_squard < ((Dneigh + Dthis) / 2 * 1.01) ** 2:
                            success = False
                            break

                shadow_pores = []
                for e in range(3):
                    if (pos[e] - Dthis/2)  < 0.0:
                        shadow_pos = np.array(pos)
                        shadow_pos[e] += bounds[e]
                        shadow_pores.append(shadow_pos)
                    if (pos[e] + Dthis/2) > bounds[e]:
                        shadow_pos = np.array(pos)
                        shadow_pos[e] -= bounds[e]
                        shadow_pores.append(shadow_pos)

                for shadow_pos in shadow_pores:
                    if tree is not  None:
                        results = tree.query_ball_point(shadow_pos, radius_min)
                        for idx in results:
                            Dneigh = diameters[idx]
                            dist_squard = ((coords[idx] - shadow_pos) ** 2).sum()
                            if dist_squard < ((Dneigh + Dthis) / 2 * 1.01) ** 2:
                                success = False
                                break
            
            if success and pos is not None:
                coords.append(pos)
                for shadow_pos in shadow_pores:
                    coords.append(shadow_pos)
                    diameters = np.insert(diameters, cidx, Dthis)
                    cidx += 1
            else:
                drops.append(cidx)
            cidx += 1

        Dmin = diameters[-1]
        Dpores = np.delete(diameters, drops)[: len(coords)]

        return PoreNetworkConfig(coords=coords, diameters=Dpores, min_diameter=Dmin)
