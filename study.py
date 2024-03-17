import numpy as np
from networks.helpers import generate_network, update_network_connections, model_from_network  # noqa
from networks.generators.config import NetworkConfig  # noqa
from networks.distribution.distribution_file import CSVFileDiameterDistribution  # noqa
from networks.generators.pores import RandomEquidistantSpacePoreGenerator, RandomPoreGenerator  # noqa
from networks.generators.conns.conns_nearest import NearestConnsGenerator
import openpnm as op
import json
from pathlib import Path
import itertools
from multiprocessing import Pool
from datetime import timedelta
import time

start = time.time()

specimens = ["O1_50", "O2_40", "O4_40"]
num_pores = [3000] # [100, 500, 1000, 1500, 2000, 2500, 3000, 3500]

ratio_throat_pores = [
    0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 
    0.01,  0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.1, 0.2, 0.3, 0.4, 0.5
]
max_tries = 10
num_variants = 50

combinations = [specimens, num_pores, ratio_throat_pores]
combinations = list(itertools.product(*combinations))

study_result_file = (Path.cwd() / "data/studies/study.json")

phys = op.models.collections.physics.basic
del phys['throat.entry_pressure']


with (Path.cwd() / f"data/data.json").open("r") as fp:
    data = json.load(fp)


def study(specimen, num_pore, ratio):
    distrib_file = './data/poresizes/comulative_pore_volume/' + specimen + '.csv'
    Deff_exact = data[specimen]["conductivity"]
    porosity = data[specimen]["porosity"]
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
            electrolyte["throat.diffusivity"] = data[data[specimen]["conductor"]]["conductivity"] #0.00023859403095188725# [S/cm] ionic conductivity
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
    print(num_pore, ":", ratio, ", ERRROR:", err)
    return err


with study_result_file.open("r+") as fp:
    results = json.load(fp)

for i in range(len(combinations)):
    combi = combinations[i]
    specimen, num_pore, ratio = combi[0], combi[1], combi[2]
    
    if results.get(specimen) is None:
        results[specimen] = {}
    
    if results[specimen].get(int(num_pore)) is None:
        results[specimen][int(num_pore)] = {}

    if results[specimen][int(num_pore)].get(str(ratio)) is None:
        results[specimen][int(num_pore)][str(ratio)] = []
    
    def process(id):
        return study(specimen, num_pore, ratio)
    
    try:
        variants = results[specimen][int(num_pore)][str(ratio)]
        with Pool() as pool: 
            result = pool.map(process, range(len(variants),num_variants))
        results[specimen][int(num_pore)][str(ratio)] = variants + result
    finally:
        with study_result_file.open("w+") as fp:
            json.dump(results, fp)
    
    elapsed = (time.time() - start)
    print(f"finished {i}/{len(combinations)}: elapsed time: " + str(timedelta(seconds=elapsed)))

print("Study ended successfully!")

