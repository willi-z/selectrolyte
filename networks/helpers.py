from .generators.config import NetworkConfig
from .distribution.idistribution import IDiameterDistribution
from .generators.pores.ipore import IPoreGenerator
from .generators.pores.helpers import get_pore_vol
from .generators.conns.iconns import IConnsGenerator
from .generators.conns.config import ConnsNetworkConfig
import numpy as np

from .config import NetworkModelConfig
import openpnm as op



def generate_network(
    diameter_generator: IDiameterDistribution,
    num_pores: int,
    porosity: float, 
    vol_ratio_throats_spheres: float,
    pore_generator: IPoreGenerator,
    conns_generator: IConnsGenerator
    ) -> NetworkConfig:
    assert 0 < vol_ratio_throats_spheres < 1.0
    pore_diameters = diameter_generator.generate_n_diameters(num_pores)
    pore_vol = get_pore_vol(pore_diameters)
    throat_vol = pore_vol * vol_ratio_throats_spheres
    box_vol = (pore_vol + throat_vol) / porosity
    box_length = np.power(box_vol,1/3.)
    bounds = (box_length, box_length, box_length)
    # print("[generate_network] vol_empty:", pore_vol + throat_vol)
    # print("[generate_network] vol_pores:", pore_vol)
    # print("[generate_network] vol_throats:", throat_vol)
    # print(vol_ratio_throats_spheres, box_vol, box_length) # <---
    pore_config = pore_generator.generate(pore_diameters, bounds)
    conns_config = ConnsNetworkConfig(conns=[], diameters=[])
    if conns_generator is not None:
        conns_config = conns_generator.generate(bounds, porosity, pore_config)
    return NetworkConfig(pores=pore_config, conns=conns_config, bounds=bounds, porosity=porosity)


def update_network_connections(
    conf: NetworkConfig, 
    conns_generator: IConnsGenerator
    )-> NetworkConfig:
    porosity = conf.porosity
    conns_config = conns_generator.generate(conf.bounds, porosity, conf.pores)
    return NetworkConfig(pores=conf.pores, conns=conns_config, bounds=conf.bounds, porosity=conf.porosity)


def model_from_network(
    conf: NetworkConfig,
    eps:float = 0.05
    )-> NetworkModelConfig:

    net = op.network.Network(coords=conf.pores.coords, conns=conf.conns.conns)
    # print(len(coords) == len(Dpores), len(coords), len(Dpores))
    net['pore.diameter'] = conf.pores.diameters
    net['throat.diameter'] = conf.conns.diameters

    net.add_model(propname='pore.volume',
                 model=op.models.geometry.pore_volume.sphere)
    net.add_model(propname='throat.length',
                 model=op.models.geometry.throat_length.spheres_and_cylinders)
    net.add_model(propname='throat.total_volume',
                 model=op.models.geometry.throat_volume.cylinder)
    net.add_model(propname='throat.lens_volume', 
                model=op.models.geometry.throat_volume.lens)
    net.add_model(propname='throat.volume', 
                 model=op.models.misc.difference,
                 props=['throat.total_volume', 'throat.lens_volume'])


    net.add_model(
        propname='throat.hydraulic_size_factors', 
        model=op.models.geometry.hydraulic_size_factors.spheres_and_cylinders
    )

    net.add_model(
        propname='throat.diffusive_size_factors', 
        model=op.models.geometry.diffusive_size_factors.spheres_and_cylinders
    )
    net.regenerate_models()
    xPore = np.array(conf.pores.coords)[:, 0]

    inner_pores = np.ones(len(conf.pores.coords))
    left_pores = np.full(len(conf.pores.coords), False, dtype=bool)
    right_pores = np.full(len(conf.pores.coords), False, dtype=bool)
    x_half = conf.bounds[0] / 2

    for i in range(len(conf.pores.coords)):
        pos = conf.pores.coords[i]
        if np.any(np.less(pos, np.zeros(3))):
            inner_pores[i] = 0
        elif np.any(np.greater(pos, conf.bounds)):
            inner_pores[i] = 0
        else:
            if pos[0] < x_half:
                left_pores[i] = True
            else:
                right_pores[i] = True
    net["pore.inner"] = inner_pores

    throat_impact = np.ones(len(conf.conns.conns))
    for i in range(len(throat_impact)):
        conns = conf.conns.conns[i]
        if inner_pores[conns[0]] == 0:
            throat_impact[i] -= 0.5
        if inner_pores[conns[1]] == 0:
            throat_impact[i] -= 0.5
    net['throat.volume_impact'] = throat_impact
    net["pore.left"] = np.isclose(xPore, 0.0, atol=conf.bounds[0]*eps)
    net["pore.right"] = np.isclose(xPore, conf.bounds[0], atol=conf.bounds[0]*eps)

    net["pore.left"] = left_pores
    net["pore.right"] = right_pores

    # print("left: ", np.where(net['pore.left'])[0])
    # print("right: ", np.where(net['pore.right'])[0])

    Vol_void = np.sum(net['pore.volume'] * inner_pores)  +np.sum(net['throat.volume'] * net['throat.volume_impact'])
    Vol_bulk = conf.bounds[0] * conf.bounds[1] * conf.bounds[2]
    porosity = Vol_void / Vol_bulk
    print(f'Porosity is: {porosity * 100:.2f}% (goal was: {conf.porosity* 100:.2f}%)')
    return NetworkModelConfig(network=net, bounds=conf.bounds, porosity=porosity)