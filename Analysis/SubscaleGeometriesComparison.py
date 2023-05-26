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
        'thrust': 1200 * const.lbToN,  # Newtons
        'time': 7.5,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 10 * const.psiToPa,
       # 'phi':1,
        'cr': 6,
        'lstar': 1.24,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 *2,
        'dp': 150 * const.psiToPa,
        'impulseguess' :  495555.24828424345,
        #'rc' : .11,
        'thetac' : (35*math.pi/180),
        'isp_efficiency' : .9} #623919}



for cr in [5,11]:
    for lstar in [1,1.24]:
        args['cr'] = cr
        args['lstar'] = lstar
        ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=False)
        args['rm']=mrideal
        params = FAC.SpreadsheetSolver(args)

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
                                        params['re'], params['lc']*.5, .0254, .0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
        # xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar
        plt.plot(np.hstack((np.flip(xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) ,label=f"lstar = {params['lstar']}, cr = {params['cr']}")

plt.legend()
plt.title(f"Different Geometries for subscale with {params['thrust']/const.lbToN} lb thrust")

plt.figure()
for pe in [10*const.psiToPa,14.7*const.psiToPa]:
    args['cr'] = 5.295390217
    args['lstar'] = 1
    args['pe'] = pe
    ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=False)
    args['rm']=mrideal
    params = FAC.SpreadsheetSolver(args)

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
                                    params['re'], params['lc']*.5, .0254, .0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
    # xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar
    plt.plot(np.hstack((np.flip(xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) ,label=f"er = {params['er']}, cr = {params['cr']}, pe = {params['pe']}")

plt.legend()
plt.title(f"Different Geometries varying exit pressure for subscale with {params['thrust']/const.lbToN} lb thrust")
plt.grid()
plt.show()