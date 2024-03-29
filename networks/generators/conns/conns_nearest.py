from .iconns import IConnsGenerator
from .config import ConnsNetworkConfig
from ..pores.config import PoreNetworkConfig
from ..pores.helpers import get_pore_vol
from scipy.spatial import KDTree
from .helpers import calc_porosity, interval_method
import numpy as np
from typing import Callable
import random


class NearestConnsGenerator(IConnsGenerator):
    def __init__(
        self, 
        num_neighbours: int = 4, # 3, 
        dfunc: Callable[[float, float, float], float] = 
        lambda d0, d1, distance: min(d0, d1)
        ):
        self.num_neighbours = num_neighbours
        self.dfunc = dfunc

    def generate(self, bounds:tuple[float, float, float], porosity: float, pores: PoreNetworkConfig) -> ConnsNetworkConfig:
        tree = KDTree(
                pores.coords
        )
        tubes = {}
        Dmax = np.max(pores.diameters)
        
        # print(vol_pores)
        box_vol = bounds[0] * bounds[1] * bounds[2]
        vol_empty = box_vol * porosity
        vol_throats = vol_empty - pores.vol
        #print("[NearestConns] vol_empty", vol_empty)
        #print("[NearestConns] vol_pores", vol_pores)
        #print("[NearestConns] vol_throats:", vol_throats)
        # assert False
        if vol_throats <= 0:
            print("vol_throats", vol_throats)
        assert vol_throats > 0

        for i in range(len(pores.coords)):
            pos = np.array(pores.coords[i])
            Dthis = pores.diameters[i]

            if True:
                if self.num_neighbours is not None:
                    distances, nidxs = tree.query(pos, self.num_neighbours)
                    if self.num_neighbours == 1:
                        nidxs = [nidxs]
                        distances = [distances]

                    for j in random.sample(range(len(nidxs)), min(len(nidxs), self.num_neighbours)):
                        idx = nidxs[j]
                        distance = distances[j]
                        if idx == i:
                            continue
                        Dneigh = pores.diameters[idx]
                        DTube = self.dfunc(Dthis, Dneigh, distance)
                        tubes[frozenset([idx, i])]=DTube
                    

            if False:
                radius_max = (pores.diameters[i] + Dmax)/2 * 1.4
                results = tree.query_ball_point(pos, radius_max)
                neighbours = set(results)
                neighbours.remove(i)
                # print(len(neighbours))
                if neighbours is None:
                    continue
                for idx in random.sample(neighbours, min(len(neighbours), self.num_neighbours)):
                    Dneigh = pores.diameters[idx]
                    pos_neigh = np.array(pores.coords[idx])
                    distance = np.sqrt(((pos_neigh - pos)**2).sum())
                    #distance_max = (Dneigh + Dthis)/2 * 1.2
                    #if distance > distance_max:
                        # continue
                    DTube = self.dfunc(Dthis, Dneigh, distance)
                    tubes[frozenset([idx, i])]=DTube

        conns = []
        diameters = []
        for conn, dia in tubes.items():
            conns.append(list(conn))
            diameters.append(dia)

        diameters = np.array(diameters)

        throat_impact = np.ones(len(conns))
        # for i in range(len(throat_impact)):
        #     conn = conns[i]
        #     if inner_pores[conn[0]] == 0.0:
        #         throat_impact[i] -= 0.5
        #     if inner_pores[conn[1]] == 0.0:
        #         throat_impact[i] -= 0.5

        def porosty_optimizer(scale):
            return calc_porosity(
                pores.coords, pores.diameters, 
                conns, diameters * scale, 
                box_vol, pores.vol,
                throat_impact) - porosity

        scale, poro = interval_method(porosty_optimizer, (0.0, 1.0), 1e-5)
        print("arch. porosity:",poro+porosity, "(vs.", porosity , ")", "by scaling with:", scale, "box vol:", box_vol)
        return ConnsNetworkConfig(conns= conns, diameters = diameters * scale, )
