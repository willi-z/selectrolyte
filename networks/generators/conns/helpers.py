import openpnm as op
import numpy as np


def calc_porosity(
    pores_coords: list[list[float]],
    pores_dia: list[float],
    conns:  list[tuple[int, int]],
    conns_dia: list[float],
    vol_bounds: float
    )-> float:
    net = op.network.Network(coords=pores_coords, conns=conns)
    # print(len(coords) == len(Dpores), len(coords), len(Dpores))
    net['pore.diameter'] = pores_dia
    net['throat.diameter'] = conns_dia

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

    Vol_void = np.sum(net['pore.volume'])+np.sum(net['throat.volume'])
    Vol_bulk = vol_bounds
    return Vol_void / Vol_bulk


def interval_method(func, bounds: tuple, eps) -> tuple[float, float]:
    x0 = bounds[0]
    x1 = bounds[1]
    xm = (x0 + x1)/2
    f0 = func(bounds[0])
    fm = func((bounds[0] + bounds[1] )/2)
    f1 = func(bounds[1])
    if x1 - x0 < eps:
        vals = [
            (abs(f0), x0),
            (abs(fm), xm),
            (abs(f1), x1),
        ]
        best_val = sorted(vals, key= lambda pair: pair[0])[0]
        return (best_val[1], best_val[0]) # (x, f)
    if fm * f0 < 0: # different sign
        return interval_method(func, (x0, xm), eps)
    else:
        return interval_method(func, (xm, x1), eps)
