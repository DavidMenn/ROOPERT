"""This is gonna be where phase 3 lies
all the different components are going to be called through this
Should be able to function as a component in two kinds of analysis:
optimization of time-static parameters (alittiude, weight, etc.)
simulation of time-dynamic parameters (instability? flight dynamics, etc. """

import scipy.optimize
#from Components.ThrustChamber import ThrustChamber
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer, add_new_propellant
import numpy as np
import os
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
import Components.ThrustChamber as ThrustChamber
import Components.StructuralApproximation as SA
from scipy.optimize import minimize_scalar

args = {
    'thrust': 5000 * const.lbToN,  # Newtons
    'time': 30,  # s
    # 'rho_ox' : 1141, #Kg/M^3
    # 'rho_fuel' : 842,
    'pc': 700 * const.psiToPa,
    'pe': 14.7 * const.psiToPa,
    'pambient': 14.7 * const.psiToPa,
    'phi': 1,
    # TYPICALLY WE SET 'phi'
    'cr': 4,
    'lstar': 1,
    'fuelname': 'Ethanol',
    'oxname': 'N2O',
    'throat_radius_curvature': .0254 / 2,
    'dp': 100 * const.psiToPa}

params = FAC.SpreadsheetSolver(args)
with open(os.path.join("..", "Configs", "rocket_output.txt"), "a") as f:
    for param in list(params):
        if params[param] is None:
            print(param + f", NONE", file=f)
        else:
            # try:
            #    print(param + f", " +'%.3f'%(params[param]), file=f)
            # except:
            print(param + f", {params[param]}", file=f)

xlist = np.linspace(0, params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], 10000)

rlist = RListGenerator.sharpRList(xlist, params['lc'], params['rc'], params['lc']+(params['rc']-params['rt'])/(2**.5),
                                  params['rt'], params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], params['re']) #xlist, xns, rc, xt, rt, xe, re
TC = ThrustChamber.ThrustChamber(rlist,xlist)
print(TC.rt, TC.xt, TC.xns)

machlist,preslist,templist = TC.flowSimple(params)

print(machlist[499]-machlist[498])

title="testing flowSimple for thrustchamber"
plt.figure(0)
fig1, axs1 = plt.subplots(4,1)

fig1.suptitle(title)

axs1[0].plot(TC.xlist,machlist , 'g')  # row=0, column=0
axs1[1].plot(TC.xlist,preslist , 'r')# row=1, column=0
axs1[2].plot(TC.xlist,templist , 'b') # row=0, column=0
axs1[3].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0


axs1[0].set_title('Mach')
axs1[1].set_title('Pressure')
axs1[2].set_title('temperature')
plt.savefig("Plots/thrust_chamber.png")
plt.show()




'''TO DO:
Get Tank volumes
- structural approx
- max accelerations
- sizes of engine components
- weights, shock force?
- apogee height
'''