from networks.network import op, net, A, L, porosity

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

import numpy as np
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
print("relative error: ", abs(D_eff - exact) / exact)
D_AB = electrolyte['pore.diffusivity'][0]
tau = porosity * D_AB / D_eff
print('The tortuosity is:', "{0:.6E}".format(tau))

