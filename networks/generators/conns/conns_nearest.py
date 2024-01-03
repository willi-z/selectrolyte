from .iconns import IConnsGenerator
from .config import ConnsNetworkConfig
from ..pores.config import PoreNetworkConfig
from ..pores.helpers import get_pore_vol
from scipy.spatial import KDTree
from .helpers import calc_porosity, interval_method
import numpy as np


class NearestConnsGenerator(IConnsGenerator):
    def __init__(self, num_neighbours: int = 3):
        self.num_neighbours = num_neighbours

    def generate(self, box_vol:float, porosity: float, pores: PoreNetworkConfig) -> ConnsNetworkConfig:
        tree = KDTree(
                pores.coords
        )
        tubes = {}
        Dmax = np.max(pores.diameters)
        vol_pores = get_pore_vol(np.array(pores.diameters))
        vol_empty = box_vol * porosity
        vol_throats = vol_empty - vol_pores
        assert vol_throats > 0

        for i in range(len(pores.coords)):
            pos = np.array(pores.coords[i])
            Dthis = pores.diameters[i]

            if self.num_neighbours is not None:
                _, results = tree.query(pos, self.num_neighbours)
                neighbours = set(results)
                neighbours.remove(i)
                
                for idx in neighbours:
                    Dneigh = pores.diameters[idx]
                    DTube = min(Dneigh, Dthis)
                    tubes[frozenset([idx, i])]=DTube
                    

            radius_max = (pores.diameters[i] + Dmax)/2 *1.6
            results = tree.query_ball_point(pos, radius_max)
            neighbours = set(results)
            neighbours.remove(i)
            # print(len(neighbours))
            if neighbours is None:
                continue
            for idx in neighbours:
                Dneigh = pores.diameters[idx]
                pos_neigh = np.array(pores.coords[idx])
                distance = np.sqrt(((pos_neigh - pos)**2).sum())
                distance_max = (Dneigh + Dthis)/2 * 1.5
                if distance > distance_max:
                    continue
                DTube = min(Dneigh, Dthis)
                tubes[frozenset([idx, i])]=DTube

        conns = []
        diameters = []
        for conn, dia in tubes.items():
            conns.append(list(conn))
            diameters.append(dia)

        diameters = np.array(diameters)

        def porosty_optimizer(scale):
            return calc_porosity(pores.coords, pores.diameters, conns, diameters * scale, box_vol)-porosity

        scale = interval_method(porosty_optimizer, (0.0, 1.0), 1e-5)
        return ConnsNetworkConfig(conns= conns, diameters = diameters * scale)
