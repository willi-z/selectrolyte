import numpy as np
from networks.helpers import generate_network, update_network_connections, model_from_network  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import openpnm as op
import json
from pathlib import Path

num_pore = 100
porosity = 0.54
ratio_throat_pores = np.arange(0.01, 0.5, 0.01)
max_tries = 10

distrib_file = './data/pore_distr_data.csv'

phys = op.models.collections.physics.basic
del phys['throat.entry_pressure']



def study_best_ratio(ratio):
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

    err = ratio_optimizer(ratio)
    print(num_pore, ":", ratio)
    return err


with (Path.cwd() / "data/studies/err_over_ratio.json").open("r+") as fp:
    results = json.load(fp)

iterations= {}
if results.get(int(num_pore)) is not None:
    iterations = results[int(num_pore)]


try:
    for ratio in ratio_throat_pores:
        variants = []
        if iterations.get(ratio) is not None:
            variants = iterations[ratio]
        for _ in range(len(variants),10):
            result = study_best_ratio(ratio)
            variants.append(result)
        iterations[ratio] = variants
finally:
    results[int(num_pore)] = iterations
    with (Path.cwd() / "data/studies/err_over_ratio.json").open("w+") as fp:
        json.dump(results, fp)

print("Study ended successfully!")

