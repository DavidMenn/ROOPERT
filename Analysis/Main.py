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
        'thrust': 3000 * const.lbToN,  # Newtons
        'time': 50,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi': 1.12,
        'mdot': 18,
        'fuelname': 'Ethanol_75',
        'oxname': 'LOX',
        'throat_radius_curvature': .02}

    params = FAC.SpreadsheetSolver(args)
    print(params)
    print(params['er'])
    dt=.05
    L, hlist, vlist, thrustlist, isplist = RE.rocketEquationCEA(params, mi = None, thrust = params['thrust'], burntime = params['time'], L = None, H = 100000, dt=dt, Af=None, ispcorrection = SIZINGCORRECTIONFACTOR)

    # Create a Figure with 2 rows and 2 columns of subplots:
    fig, ax = plt.subplots(2, 2)
    print(L)
    x = np.linspace(0, 5, 100)

    # Index 4 Axes arrays in 4 subplots within 1 Figure:
    ax[0, 0].plot(np.arange(0, dt * hlist.size, dt), hlist, 'g')  # row=0, column=0
    ax[1, 0].plot(np.arange(0, dt * hlist.size, dt), vlist, 'b')  # row=1, column=0
    ax[0, 1].plot(np.arange(0,  params['time'], dt), thrustlist[0:np.size(np.arange(0, params['time'], dt))], 'r')  # row=0, column=1
    ax[1, 1].plot(np.arange(0, params['time'], dt), isplist[0:np.size(np.arange(0, params['time'], dt))], 'k')  # row=1, column=1

    ax[0,0].set_title('Height')
    ax[1, 0].set_title('Velo (M/s)')
    ax[0, 1].set_title('Thrust (Newtons)')
    ax[1, 1].set_title('Isp (sec)')
    plt.show()

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