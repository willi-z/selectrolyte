from .ipore import IPoreGenerator, PoreNetworkConfig
from scipy.spatial import KDTree
import numpy as np
from .helpers import iso_2_coord, get_pore_vol


class RandomPoreGenerator(IPoreGenerator):
    def __init__(self) -> None:
        pass

    def generate(
        self, diameters: list[list[float]], bounds: tuple[float, float, float]
    ) -> PoreNetworkConfig:
        diameters[::-1].sort()
        volume = 0.0

        
        Dmax = diameters[0]
        Dmin = diameters[-1]

        drops = []
        # np.random.seed(0)
        coords = [iso_2_coord(np.random.rand(3), ((0, 0, 0), bounds))]
 
        # grid = VoxelGrid(bounds, Dmin)
        # grid.remove_filled_voxels(coords[0], diameters[0])
        cidx = 1

        while cidx < len(diameters):
            tree = KDTree(
                coords,
            )
            Dthis = diameters[cidx]
            radius_min = (Dmax + Dthis) / 2 * 1.2
            # radius = box_length
            success = False
            pos = None
            shadow_pos = None # in case an overhang occures

            counter = 0

            while counter < 1000 and not success:
                counter += 1
                success = True
                iso_coord = np.random.rand(3)
                pos = iso_2_coord(iso_coord, ((0, 0, 0), bounds))

                results = tree.query_ball_point(pos, radius_min)
                for idx in results:
                    Dneigh = diameters[idx]
                    dist_squard = ((coords[idx] - pos) ** 2).sum()
                    if dist_squard < ((Dneigh + Dthis) / 2 * 1.1) ** 2:
                        success = False
                        break

                shadow_pos = None
                has_overhang = False
                offset = np.zeros(3)
                overhang_lower = np.array(iso_coord) - Dthis
                axis_wise_overhang_lower = np.less(overhang_lower, np.zeros(3))
                if np.any(axis_wise_overhang_lower):
                    has_overhang = True
                    for e in range(3):
                        if axis_wise_overhang_lower[e]:
                            offset[e] = - bounds[e]

                overhang_upper = np.array(bounds) - (np.array(iso_coord) + Dthis)
                axis_wise_overhang_upper = np.less(overhang_upper, np.zeros(3))
                if np.any(axis_wise_overhang_upper):
                    has_overhang = True
                    for e in range(3):
                        offset[e] = + bounds[e]

                if has_overhang:
                    shadow_pos = np.array(iso_coord) + np.array(offset)
                    pos = iso_2_coord(shadow_pos, ((0, 0, 0), bounds))

                    results = tree.query_ball_point(pos, radius_min)
                    for idx in results:
                        Dneigh = diameters[idx]
                        dist_squard = ((coords[idx] - pos) ** 2).sum()
                        if dist_squard < ((Dneigh + Dthis) / 2 * 1.1) ** 2:
                            success = False
                            break
            
            if success and pos is not None:
                coords.append(pos)
                if shadow_pos is not None:
                    coords.append(shadow_pos)
                    np.insert(diameters, cidx, Dthis)
                    volume += get_pore_vol(np.array([diameters[cidx]]))
                    cidx += 1
            else:
                drops.append(cidx)
            cidx += 1

        Dmin = diameters[-1]
        Dpores = np.delete(diameters, drops)[: len(coords)]

        return PoreNetworkConfig(coords=coords, diameters=Dpores, min_diameter=Dmin, extra_volume=volume)
