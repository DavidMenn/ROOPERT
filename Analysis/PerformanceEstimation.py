"""Used to determine rocket parameters to get to target height given mass approx
Currently doesnt have a clean output, so in order to actually get the weights
of individual components and parameters you have to enter the debugger, put
a breakpoint at the structuraapprox call after the rocketequation call
and step through structural approx until it gives you all the weights
REALLY NEEDS A CLEAN OUTPUT SOON!"""
import sys
sys.path.insert(1,"./")
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer, add_new_propellant
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
import Toolbox.RListGenerator as RListGenerator
import Toolbox.RocketCEAAssister as RA
import Toolbox.IsentropicEquations as IE
import Toolbox.RocketEquation as RE
import difflib
import re as regex
from rocketprops.rocket_prop import get_prop
import Toolbox.Constant as const
import matplotlib.pyplot as plt
import FirstOrderCalcs as FAC
import Components.StructuralApproximation as SA
import json
import os
def SizeRocket_MassApprox():
    SIZINGCORRECTIONFACTOR = .9
    """ args = {
        'thrust': 5000 * const.lbToN,  # Newtons
        'time': 30,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi': 1,
        'cr' :4,
        'fuelname': 'Ethanol',
        'oxname': 'LOX',
        'throat_radius_curvature': .02}"""
    args = {
        'thrust': 1200 * const.lbToN,  # Newtons
        'time': 7.5,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'pambient' : 14.7*const.psiToPa,
        'phi': 1,
        # TYPICALLY WE SET 'phi'
        'cr': 4,
        'lstar': 1,
        'fuelname': 'Ethanol',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 100 * const.psiToPa}

    params = FAC.SpreadsheetSolver(args)
    dp=150*const.psiToPa
    miinit, lambdainit, totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank = SA.mass_approx(params['pc'],dp, 12, params['rho_fuel'],params['rho_ox'], params['thrust'], params['isp'], params['time'], params['rm'])
    #os.mkdir("Configs")
    #with open(os.path.join("..", "Configs","params.json"),"w") as file:
    #    json.dump({name: str(params[name]) for name in params}, file)
    #print(params)
    thrusttoweight_approx = 4

    #json.load()
    dt=.25

    params['thrust'] = totalmass*thrusttoweight_approx*9.81 #recompute with a resonable thrust
    params = FAC.SpreadsheetSolver(params)
    miinit, lambdainit, totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank = SA.mass_approx(params['pc'], dp, 12, params['rho_fuel'], params['rho_ox'],
                                                   params['thrust'], params['isp'], params['time'], params['rm'])
    L, mi, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA_MassAprox(params, 15000*const.lbToN,4, H=100000, dp = dp,
                                                                          dt=dt, Af=None,
                                                                          ispcorrection=.95)
    SA.mass_approx(params['pc'], dp, 12, params['rho_fuel'], params['rho_ox'],
                                                   params['thrust'], params['isp'], params['time'], params['rm'])
    with open(os.path.join("Configs", "mass_output.txt"), "a") as f:
        for param in list(params):
            if params[param] is None:
                print(param + f", NONE", file=f)
            else:
                #try:
                #    print(param + f", " +'%.3f'%(params[param]), file=f)
                #except:
                print(param + f", {params[param]}", file=f)

    title = f"Fuel = {params['fuelname']}, Thrust = {params['thrust']}, burntime = {params['time']}, " \
            f"mi = {((1 / L) * params['mdot'] * params['time'] - params['mdot'] * params['time'])}, mp = {params['mdot'] * params['time']}, Mtotal = {(1 / L) * params['mdot'] * params['time']}"
    RE.ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time=params['time'], title=title, dt=dt)
    plt.show()

SizeRocket_MassApprox()
