from .generators.config import NetworkConfig
from .distribution.idistribution import IDiameterDistribution
from .generators.pores.ipore import IPoreGenerator
from .generators.pores.helpers import get_pore_vol
from .generators.conns.iconns import IConnsGenerator
import numpy as np

from .config import NetworkModelConfig
import openpnm as op



def generate_network(
    diameter_generator: IDiameterDistribution,
    num_pores: int,
    porosity: float, 
    pore_generator: IPoreGenerator,
    conns_generator: IConnsGenerator
    ) -> NetworkConfig:
    pore_diameters = diameter_generator.generate_n_diameters(num_pores)
    pore_vol = get_pore_vol(pore_diameters)

    box_vol = pore_vol / porosity
    box_length = np.power(box_vol,1/3.)
    bounds = (box_length, box_length, box_length)
    pore_config = pore_generator.generate(pore_diameters, bounds)
    conns_config = conns_generator.generate(box_vol, porosity, pore_config)
    return NetworkConfig(pores=pore_config, conns=conns_config, bounds=bounds, porosity=porosity)


def update_network_connections(
    conf: NetworkConfig, 
    conns_generator: IConnsGenerator
    )-> NetworkConfig:
    box_vol = conf.bounds[0] * conf.bounds[1] * conf.bounds[2]
    porosity = conf.porosity
    conns_config = conns_generator.generate(box_vol, porosity, conf.pores)
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
    net["pore.left"] = np.isclose(xPore, 0.0, atol=conf.bounds[0]*eps)
    net["pore.right"] = np.isclose(xPore, conf.bounds[0], atol=conf.bounds[0]*eps)

    # print("left: ", np.where(net['pore.left'])[0])
    # print("right: ", np.where(net['pore.right'])[0])

    Vol_void = np.sum(net['pore.volume'])+np.sum(net['throat.volume'])
    Vol_bulk = conf.bounds[0] * conf.bounds[1] * conf.bounds[2]
    porosity = Vol_void / Vol_bulk
    print(f'Porosity is: {porosity * 100:.2f}% (goal was: {conf.porosity* 100:.2f}%)')
    return NetworkModelConfig(network=net, bounds=conf.bounds, porosity=porosity)