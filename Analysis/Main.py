#from Components.ThrustChamber import ThrustChamber
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
#xlist = np.linspace(0, .35, 1000)
#rlist = RListGenerator.sharpRList(xlist, .2, .07, .24, .02, .35, .06) #xlist, xns, rc, xt, rt, xe, re
#TC = ThrustChamber(rlist,xlist)
#print(TC.rt, TC.xt, TC.xns)

#RA.makeEthanolBlend(75)

#C = CEA_Obj(oxName="LOX", fuelName="Ethanol_75")

#s = C.get_full_cea_output( Pc=350, MR=[1.19,1.2,1.21] ,eps=4, short_output=1)

#print( s )
def main():
    SIZINGCORRECTIONFACTOR = .9
    args = {
        'thrust': 5000 * const.lbToN,  # Newtons
        'time': 30,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi': 1,
        'fuelname': 'Ethanol',
        'oxname': 'LOX',
        'throat_radius_curvature': .02}

    params = FAC.SpreadsheetSolver(args)
    dp=150*const.psiToPa
    miinit,lambdainit,totalmass = SA.mass_approx(params['pc'],dp, 12, params['rho_fuel'],params['rho_ox'], params['thrust'], params['isp'], params['time'], params['rm'])

    print(params)
    thrusttoweight_approx = 10

    dt=.25

    params['thrust'] = totalmass*thrusttoweight_approx*9.81 #recompute with a resonable thrust
    params = FAC.SpreadsheetSolver(params)
    miinit, lambdainit, totalmass = SA.mass_approx(params['pc'], dp, 12, params['rho_fuel'], params['rho_ox'],
                                                   params['thrust'], params['isp'], params['time'], params['rm'])
    L, mi, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA_MassAprox(params, 15000*const.lbToN,4, H=100000, dp = dp,
                                                                          dt=dt, Af=None,
                                                                          ispcorrection=.95)
    SA.mass_approx(params['pc'], dp, 12, params['rho_fuel'], params['rho_ox'],
                                                   params['thrust'], params['isp'], params['time'], params['rm'])
    title = f"Fuel = {params['fuelname']}, Thrust = {params['thrust']}, burntime = {params['time']}, " \
            f"mi = {((1 / L) * params['mdot'] * params['time'] - params['mdot'] * params['time'])}, mp = {params['mdot'] * params['time']}, Mtotal = {(1 / L) * params['mdot'] * params['time']}"
    RE.ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time=params['time'], title=title, dt=dt)
    """

    print(params)
    print(params['er'])
    print(params['temp_throat'])
    dt=.05
    L, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA(params, mi = None, thrust = params['thrust'], burntime = params['time'], L = None, H = 100000, dt=dt, Af=None, ispcorrection = SIZINGCORRECTIONFACTOR)

    # Create a Figure with 2 rows and 2 columns of subplots:

    print("mi = " + str((1 / L) * params['mdot']*params['time']-params['mdot']*params['time']))
    print('mp = ' + str(params['mdot']*params['time']))
    print("Mtotal = "+ str((1 / L) * params['mdot']*params['time']))
    title = f"Fuel = {params['fuelname']}, Thrust = {params['thrust']}, burntime = {params['time']}, " \
            f"mi = {((1 / L) * params['mdot']*params['time']-params['mdot']*params['time'])}, mp = {params['mdot']*params['time']}, Mtotal = {(1 / L) * params['mdot']*params['time']}"

    RE.ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time=params['time'], title=title, dt=None)
    """
main()

"""    
    TransportCEA = params['CEA'].get_Chamber_Transport(Pc=params['pc'], MR=params['rm'], eps=params['er'],frozen=0)
    print(TransportCEA)

    #print(4*params['gamma_throat']/(9*params['gamma_throat']-5))
    #Bartz(params['rt']*2,)


def Bartz(D, viscosityns, prandtlns, cpns, pcns, cstar, throat_radius_curvature, at, a, twg, tcns, mach, gamma):
    sigma = ((.5 * twg / tcns * (1 + (mach ** 2) * (gamma - 1) / 2) + .5) ** (-.68)) * (
                1 + (mach ** 2) * (gamma - 1) / 2) ** (-.12)
    return (.026 / (D ** .2) * (cpns * viscosityns ** .2) / (prandtlns ** .6) * (pcns / cstar) ** (.8) * (
                D / throat_radius_curvature) ** .1) * (at / a) ** .9 * sigma
main()
"""