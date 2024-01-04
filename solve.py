from networks.generators.config import NetworkConfig
from networks.helpers import model_from_network
from pathlib import Path
import json
import openpnm as op
import numpy as np

with (Path.cwd() / "data/networks/rand_0.54.json").open("r") as fp:
    content = json.load(fp)

conf = NetworkConfig(**content)
model = model_from_network(conf)
net = model.network
L = model.bounds[0]
A = model.bounds[1] * model.bounds[2]

# phase
# liquid = op.phase.Water(network=net)
electrolyte = op.phase.Phase(network=net)
# electrolyte = op.phase.Air(network=net)

phys = op.models.collections.physics.basic
del phys['throat.entry_pressure']
electrolyte.add_model_collection(phys)
electrolyte["pore.viscosity"] = 0.05  # [Pa.s] seems to have no impact -> unrelevent
electrolyte["throat.diffusivity"] = 0.00023859403095188725# [S/cm] ionic conductivity
electrolyte.regenerate_models()

fd = op.algorithms.FickianDiffusion(network=net, phase=electrolyte)


if not np.isfinite(fd.A.data).all():
    print(fd.A.data)
    print("A")
if not np.isfinite(fd.b).all():
    print("b")

inlet = net.pores('left')
outlet = net.pores('right')
C_in, C_out = [10, 5]
fd.set_value_BC(pores=inlet, values=C_in)
fd.set_value_BC(pores=outlet, values=C_out)

fd.run()

rate_inlet = fd.rate(pores=inlet)[0]
print(f'Molar flow rate: {rate_inlet:.5e} mol/s')


D_eff = rate_inlet * L / (A * (C_in - C_out))
print("The effective  diffusivity is: {0:.6E} m/s^2".format(D_eff))
# 
exact = 0.000005321864708362311
print("abs. error", D_eff - exact)
print("rel. error: ", abs(D_eff - exact) / exact)
D_AB = electrolyte['pore.diffusivity'][0]
tau = model.porosity * D_AB / D_eff
print('The tortuosity is:', "{0:.6E}".format(tau))

