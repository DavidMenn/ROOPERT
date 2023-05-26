"""Just a place to run stuff, currently has some test code
usefull as it has a bunch of commonly used import statements
CURRENT CONFIGURATION: FIRST ORDER SIZINGS FOR MID AUTUMN 2022"""
import sys
sys.path.insert(1,"./")
import scipy.optimize
#from Components.ThrustChamber import ThrustChamber
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer, add_new_propellant
import numpy as np
import math
import Components.ThrustChamber as ThrustChamber
import Components.CoolingSystem as CS
from rocketcea.cea_obj_w_units import CEA_Obj
import Toolbox.RListGenerator as RListGenerator
import Toolbox.RocketCEAAssister as RA
import os
import Toolbox.IsentropicEquations as IE
import Toolbox.RocketEquation as RE
import difflib
import re as regex
from rocketprops.rocket_prop import get_prop
import Toolbox.Constant as const
import DetermineOptimalMR as DOMR
import matplotlib.pyplot as plt
import FirstOrderCalcs as FAC
import Components.ThrustChamber as ThrustChamber
import Components.StructuralApproximation as SA
from scipy.optimize import minimize_scalar
from Toolbox import PressureDropCalculator as PD


args = {
        'thrust': 5000 * const.lbToN,  # Newtons
        'time': 30,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi': 1,
        'cr' :4,
        'lstar' : 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .02}

ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=False)
print(f"isp max = {ispmaxavg}, ideal mr is {mrideal}")
args['rm']=mrideal
params = FAC.SpreadsheetSolver(args)

# get pressure drop

try:
    try:
        params['dp'] = [params['fueldp'],params['oxdp']]
    except:
        params['dp'] = [params['dp'],params['dp']]
except: #mess of a try except statement to find pressure drops, if you specify dp itll just do that for both. if you don't, itll try to use the pd calculator
    try:
        oxpd, fuelpd, dparray = PD.pressureDrop(OxDensityMetric = params['rho_ox'],     # kg/m^3
                                    OxMassFlowMetric = params['mdot_ox'], # kg/s
                                    OxTubeDiam = params['OxTubeDiam'], #in  
                                    FuelDensityMetric = params['rho_fuel'],       # kg/m^3
                                    FuelMassFlowMetric = params['mdot_fuel'],# kg/s
                                    FuelTubeDiam = params['FuelTubeDiam'],
                                    params = params )
        try: 
            params['dp'] =  [(params['pc']*params['injectorstiffness'] + params['regenpd'] + fuelpd)*params['pdsafetyfactor'], (params['pc']*params['injectorstiffness'] + oxpd)*params['pdsafetyfactor']]
        except:
            params['dp'] =  [(params['pc']*.2 + 50*const.psiToPa + fuelpd)*1.5, (params['pc']*.2 + oxpd)*1.5]
    except:
        params['dp'] =  [150*const.psiToPa ,150*const.psiToPa ]


if params['thetac'] is None:
    params['thetac'] = math.pi*35/180
volfunc = lambda lc : math.pi*params['rc']**2*lc  +\
    math.pi*params['rc']**3/math.tan(params['thetac'])/3 -\
        math.pi*params['rt']**3/math.tan(params['thetac'])/3
lstarminimizer = lambda lc : volfunc(lc)/(params['rt']**2*math.pi) - params['lstar']
result = scipy.optimize.root(lstarminimizer, .05, args=(), method='hybr', jac=None, tol=None, callback=None, options=None)
params['lc']=result['x'][0]
xlist = np.linspace(0, params['lc'] + (params['rc'] - params['rt']) / math.tan(params['thetac']) + params['ln_conical'], 100)
    
rlist,xlist = RListGenerator.paraRlist(xlist, params['lc'], params['rc'],
                                params['lc'] + (params['rc'] - params['rt'])/(math.tan(params['thetac'])),
                                params['rt'],
                                params['lc'] + (params['rc'] - params['rt'])/(math.tan(params['thetac'])) + params['ln_conical'],
                                params['re'], params['lc']*.5, 2*.0254, 2*.0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
# xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar

TC = ThrustChamber.ThrustChamber(rlist,xlist)
print(TC.rt, TC.xt, TC.xns)

machlist,preslist,templist = TC.flowSimple(params)


title="FlowSimple Test"
fig, axs = plt.subplots(2,2)
fig.suptitle(title)

axs[0,0].plot(xlist,templist, 'g')  
axs[0,1].plot(xlist,preslist, 'r')
axs[1,0].plot(xlist,machlist, 'b') 
axs[1,1].plot(xlist,math.pi*rlist**2)

axs[0,0].set_title('templist')
axs[0,1].set_title('preslist')
axs[1,0].set_title('machlist')
axs[1,1].set_title('arealist')
plt.show()

#machtest = IE.machFromArea(a,at,gam,supersonic = False):
print(IE.machFromArea(2,1,1.4,supersonic = False)) #A/At = 2
print(IE.machFromArea(3,1,1.4,supersonic = True)) #A/At = 3
print(IE.machFromArea(4,1,1.4,supersonic = False)) #A/At = 4
print(IE.machFromArea(5,1,1.4,supersonic = True)) #A/At = 5
