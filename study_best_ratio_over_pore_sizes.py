import numpy as np
from networks.helpers import generate_network, update_network_connections, model_from_network  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import openpnm as op
import json
from pathlib import Path

num_pores = np.arange(100, 200, 2)
porosity = 0.54
ratio_throat_pores = 0.02
max_tries = 10

distrib_file = './data/pore_distr_data.csv'

phys = op.models.collections.physics.basic
del phys['throat.entry_pressure']


def interval_method(func, bounds: tuple[float, float], eps) -> tuple[float, float]:
    x0 = bounds[0]
    x1 = bounds[1]
    xm = (x0 + x1)/2
    f0 = func(bounds[0])
    fm = func((bounds[0] + bounds[1] )/2)
    f1 = func(bounds[1])
    if f1 * f0 > 0:
        print("ERROR",f0, "vs.", f1, "(should be differnt signs)")
        print("Arguments",x0, "and", x1)
        assert f1 * f0 <= 0 # different sign
    if x1 - x0 < eps:
        vals = [
            (x0, f0),
            (xm, fm),
            (x1, f1),
        ]
        best_val = sorted(vals, key= lambda pair: abs(pair[1]))[0]
        return (best_val[0], best_val[1]) # (x, f)
    if fm * f0 <= 0: # different sign
        return interval_method(func, (x0, xm), eps)
    else:
        return interval_method(func, (xm, x1), eps)


def study_best_ratio(num_pore):
    Deff_exact = 0.000005321864708362311
    def ratio_optimizer(ratio):
        num_tries = 0
        success = False
        while num_tries < max_tries and not success:
            conf = generate_network(
                CSVFileDiameterDistribution(distrib_file),
                num_pore,
                porosity,
                ratio,
                RandomPoreGenerator(),
                NearestConnsGenerator()
            )
            model = model_from_network(conf)
            net = model.network
            L = model.bounds[0]
            A = model.bounds[1] * model.bounds[2]

            # phase
            # liquid = op.phase.Water(network=net)
            electrolyte = op.phase.Phase(network=net)
            # electrolyte = op.phase.Air(network=net)


            electrolyte.add_model_collection(phys)
            electrolyte["pore.viscosity"] = 0.05  # [Pa.s] seems to have no impact -> unrelevent
            electrolyte["throat.diffusivity"] = 0.00023859403095188725# [S/cm] ionic conductivity
            electrolyte.regenerate_models()

            try:
                fd = op.algorithms.FickianDiffusion(network=net, phase=electrolyte)

                inlet = net.pores('left')
                outlet = net.pores('right')
                C_in, C_out = [10, 5]
                fd.set_value_BC(pores=inlet, values=C_in)
                fd.set_value_BC(pores=outlet, values=C_out)

                fd.run()
            except Exception:
                num_tries += 1
                continue
            success = True

        if not success:
            raise ValueError

        rate_inlet = fd.rate(pores=inlet)[0]

        D_eff = rate_inlet * L / (A * (C_in - C_out))

        
        return D_eff - Deff_exact

    eps = 1e-4
    ratio, err = interval_method(ratio_optimizer, (eps, 0.1), eps)
    print(num_pore, ":", ratio)
    return dict(ratio=ratio, err=err)


with (Path.cwd() / "data/studies/best_ratio_over_poresizes.json").open("r+") as fp:
    results = json.load(fp)

try:
    for num_pore in num_pores:
        if results.get(num_pore) is not None:
            continue
        result = study_best_ratio(num_pore)
        results[int(num_pore)] = result
finally:
    with (Path.cwd() / "data/studies/best_ratio_over_poresizes.json").open("w+") as fp:
        json.dump(results, fp)

print("Study ended successfully!")

