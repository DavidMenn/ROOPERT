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

def fixedRSquareChanelSetup(params,xlist, rlist,chlist,chanelToLandRatio,twlist,nlist,helicitylist = None,dxlist = None):# MAKE SURE TO PASS SHIT  ALREADY FLIPPED
    if helicitylist is None:
        helicitylist = np.ones(np.size(xlist))*math.pi/2
    if dxlist is None:
        dxlist=np.ones(np.size(xlist))
        index = 1
        while index<np.size(xlist):
            dxlist[index]=abs(xlist[index-1]-xlist[index])
            index = index+1
        dxlist[0]=dxlist[1] # this is shit, but I actuall calculate an extra spatial step at the start for some reason, so our CC is 1 dx too long. Makes the graphs easier to work with tho lol, off by one error be damned

    #FOR NOW THIS IS ASSUMING SQUARE CHANELS!
    #\     \  pavail\     \ 
    # \     \ |----| \     \
    #  \     \        \     \
    #   \     \        \     \
    # sin(helicitylist) pretty much makes sure that we are using the paralax angle to define the landwidth
    perimavailable = math.pi*2*(rlist+twlist)/nlist*np.sin(helicitylist)
    landwidthlist=perimavailable/(chanelToLandRatio+1)
    cwlist=perimavailable-landwidthlist

    alistflipped=chlist*cwlist
    salistflipped=cwlist/dxlist
    vlistflipped = params['mdot_fuel'] / params['rho_fuel'] / alistflipped / nlist
    #if np.min(cwlist/chlist)>10:
    #    raise Exception(f"Aspect Ratio is crazyyyyy")

    hydraulicdiamlist=4*alistflipped/(2*chlist+2*cwlist)
    coolingfactorlist = np.ones(xlist.size)
    heatingfactorlist = np.ones(xlist.size)*.6 # .6 is from cfd last year, i think its bs but whatever
    """Fin cooling factor func is 2*nf*CH+CW. nf is calculated as tanh(mL)/ml. Ml is calculated as sqrt(hpL^2/ka),
    h=hc, P=perimeter in contact with coolant = dx, A = dx*landwidth/2 (assuming only half the fin works on the coolant, 2* factor in other spot,
    I think I can do that because of axisymmetric type), k=kw, L=height of fin = chanel height
    This is all from "A Heat Transfer Textbook, around page 166 and onwards."""
    fincoolingfactorfunc = lambda hc,kw,ind : (math.tanh(chlist[ind]*math.sqrt(2*hc/kw*landwidthlist[ind]))/\
            (chlist[ind]*math.sqrt(2*hc/kw*landwidthlist[ind])))*2*chlist[ind] + cwlist[ind]

    return alistflipped, nlist, coolingfactorlist, heatingfactorlist, xlist, vlistflipped, twlist, hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist

def RunCoolingSystem(chlist,
twlist,
nlist,
helicitylist,
params,
xlist,
rlist,
chanelToLandRatio, 
TC ):

    machlist,preslist,templist = TC.flowSimple(params)
    alistflipped, nlist, coolingfactorlist, heatingfactorlist, xlist, vlistflipped, twlistflipped,\
     hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist   = fixedRSquareChanelSetup(params = params,
                                        xlist = np.flip(xlist), rlist = np.flip(rlist),chlist = chlist ,
                                        chanelToLandRatio = chanelToLandRatio ,twlist = twlist ,nlist = nlist,
                                        helicitylist=helicitylist)

    Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist = CS.steadyStateTemperatures(
        None, TC, params, salistflipped, nlist, coolingfactorlist,
        heatingfactorlist, xlist, vlistflipped, 293, params['pc'] + params['pc']*.2 + 50*const.psiToPa, twlistflipped, hydraulicdiamlist, rgaslist = rlist, fincoolingfactorfunc=fincoolingfactorfunc, dxlist = dxlist)

    material = "inconel 715"
    Structure = CS.StructuralAnalysis(rlist, xlist, nlist, chlist, cwlist, twlist, material)
    FOSlist = Structure.FOS(Twglist,Twclist,coolantpressurelist,preslist)
    xlist=np.flip(xlist) #idk whats flipping this haha but its something in the steadystatetemps function, so we have to flip it back
    return alistflipped, xlist, vlistflipped,\
     hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist,\
        Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist,\
            FOSlist
def FirstOrderSolver(args):
    output=False
    try: # only computs optimal inputted rm is none
        if args['rm'] == None:
            ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=output)
            args['rm']=mrideal
    except:
        ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=output)
        args['rm']=mrideal

    
    params = FAC.SpreadsheetSolver(args)

    try:    
        try:
            params['dp'] = [params['fueldp'],params['oxdp']]
        except:
            params['dp'] = [params['p'],params['dp']]
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

    params['vol_ox'] = params['mdot_ox']*params['time']/params['rho_ox']
    params['vol_fuel'] = params['mdot_fuel']*params['time']/params['rho_fuel']
    params['mass_ox'] = params['mdot_ox']*params['time']
    params['mass_fuel'] = params['mdot_fuel']*params['time']
    params['P_tank_fuel'] = params['pc'] + params['dp'][0]
    params['P_tank_ox']= params['pc'] + params['dp'][1]
    #params['TWR'] = params['thrust']/params['totalmass']
    return params
def RocketEquation_MassFixed(params):
    dt=.025
    L, hlist, vlist, thrustlist, isplist, rocketmachlist = RE.rocketEquationCEA(params, mi = params['drymass'],
                    thrust = params['thrust'], burntime = params['time'],
                    L = None, H = None, dt=dt, Af=None, 
                    ispcorrection = params['isp_efficiency'])
    
    burntimelist = np.arange(0, params['time'], dt)
    timelist = np.arange(0, np.where(hlist==0)[0][1]*dt-dt/2, dt)
    alist = np.diff(vlist) / dt/9.81 # in gs
    lists = {}
    lists['hlist'] = hlist[0:len(timelist)]
    lists['vlist'] = vlist[0:len(timelist)]
    lists['thrustlist'] = thrustlist[0:len(burntimelist)]
    lists['isplist'] = isplist[0:len(burntimelist)]
    lists['rocketmachlist'] = rocketmachlist[0:len(timelist)]
    lists['alist'] = alist[0:len(timelist)]
    lists['timelist'] = timelist
    lists['burntimelist'] = burntimelist
    return lists

def ChanelSolver(params,chanelgeometry):
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
                                    params['re'], params['lc']*1.8, 2*.0254, 2*.0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
    # xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar

    TC = ThrustChamber.ThrustChamber(rlist,xlist)
    print(TC.rt, TC.xt, TC.xns)

    machlist,preslist,templist = TC.flowSimple(params)
    xlistflipped = np.flip(xlist)
    rlistflipped  = np.flip(rlist)
    chlist  = chanelgeometry['chlist']# (TC.rt/rlistflipped)**.5*.003 
    twlist  = chanelgeometry['twlist']#(rlistflipped/TC.rt)*.001 
    nlist  = chanelgeometry['nlist']#np.ones(len(xlist))*80
    ewlist  = chanelgeometry['ewlist']#np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylist  = chanelgeometry['helicitylist']#(rlistflipped**1.5/TC.rt**1.5)*45*math.pi/180
    chanelToLandRatio = chanelgeometry['c2l']#
    for index in range(0,np.size(helicitylist )):
        if helicitylist [index]>math.pi/2:
            helicitylist [index] = math.pi/2
    alistflipped, xlist, vlistflipped,\
        hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist,\
            Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist, FOSlist = RunCoolingSystem(chlist,
    twlist,
    nlist,
    helicitylist,
    params,
    xlist,
    rlist,
    chanelToLandRatio, 
    TC )
    lists = {}
    lists['xlist'] = xlist
    lists['rlist'] = rlist
    lists['machlist'] = machlist
    lists['preslist'] = preslist
    lists['templist'] = templist
    lists['chlist']  = chlist# (TC.rt/rlistflipped)**.5*.003 
    lists['twlist']  = twlist#(rlistflipped/TC.rt)*.001 
    lists['nlist']  = nlist#np.ones(len(xlist))*80
    lists['ewlist']  = ewlist#np.ones(len(xlist))*.005
    #HELIC'ITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    lists['helicitylist']  = helicitylist
    lists['alistflipped'] = alistflipped
    lists['vlistflipped'] = vlistflipped
    lists['hydraulicdiamlist'] = hydraulicdiamlist
    lists['salistflipped'] = salistflipped
    lists['dxlist'] = dxlist
    lists['cwlist'] = cwlist
    lists['Twglist'] = Twglist
    lists['hglist'] = hglist
    lists['qdotlist'] = qdotlist
    lists['Twclist'] = Twclist
    lists['hclist'] = hclist
    lists['Tclist'] = Tclist
    lists['coolantpressurelist'] = coolantpressurelist
    lists['fincoolingfactorlist'] = fincoolingfactorlist
    lists['rholist']=rholist
    lists['viscositylist']=viscositylist
    lists['Relist']=Relist
    lists['FOSlist']=FOSlist
def Everything_MassFixed(args, chanelgeometry):

    output=False
    try: # only computs optimal inputted rm is none
        if args['rm'] == None:
            ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=output)
            args['rm']=mrideal
    except:
        pass # hopefully it calcls isp and phi in FAC

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
    params = FAC.SpreadsheetSolver(args)
    
    params['vol_ox'] = params['mdot_ox']*params['time']/params['rho_ox']
    params['vol_fuel'] = params['mdot_fuel']*params['time']/params['rho_fuel']
    params['mass_ox'] = params['mdot_ox']*params['time']
    params['mass_fuel'] = params['mdot_fuel']*params['time']
    params['P_tank_fuel'] = params['pc'] + params['dp'][0]
    params['P_tank_ox']= params['pc'] + params['dp'][1]
    
    params['totalmass'] = params['drymass']+params['vol_fuel']*params['rho_fuel']+params['vol_ox']*params['rho_ox']

    dt=.025
    params['TWR'] = params['thrust']/params['totalmass']
    L, hlist, vlist, thrustlist, isplist, rocketmachlist = RE.rocketEquationCEA(params, mi = params['drymass'],
                    thrust = params['thrust'], burntime = params['time'],
                    L = None, H = params['target'], dt=dt, Af=None, 
                    ispcorrection = params['isp_efficiency'])
    
    burntimelist = np.arange(0, params['time'], dt)
    timelist = np.arange(0, np.where(hlist==0)[0][1]*dt-dt/2, dt)
    alist = np.diff(vlist) / dt/9.81 # in gs
    
    #conevol = math.pi*params['rc']**3*math.tan(params['thetac'])/3 - math.pi*params['rt']**3*math.tan(params['thetac'])/3
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
                                    params['re'], params['lc']*1.8, 2*.0254, 2*.0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
    # xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar

    TC = ThrustChamber.ThrustChamber(rlist,xlist)
    print(TC.rt, TC.xt, TC.xns)

    machlist,preslist,templist = TC.flowSimple(params)
    xlistflipped = np.flip(xlist)
    rlistflipped  = np.flip(rlist)
    chlist  = chanelgeometry['chlist']# (TC.rt/rlistflipped)**.5*.003 
    twlist  = chanelgeometry['twlist']#(rlistflipped/TC.rt)*.001 
    nlist  = chanelgeometry['nlist']#np.ones(len(xlist))*80
    ewlist  = chanelgeometry['ewlist']#np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylist  = chanelgeometry['helicitylist']#(rlistflipped**1.5/TC.rt**1.5)*45*math.pi/180
    chanelToLandRatio = chanelgeometry['c2l']#
    for index in range(0,np.size(helicitylist )):
        if helicitylist [index]>math.pi/2:
            helicitylist [index] = math.pi/2
    alistflipped, xlist, vlistflipped,\
        hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist,\
            Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist, FOSlist = RunCoolingSystem(chlist,
    twlist,
    nlist,
    helicitylist,
    params,
    xlist,
    rlist,
    chanelToLandRatio, 
    TC )
    
    lists = {}
    lists['hlist'] = hlist
    lists['vlist'] = vlist
    lists['thrustlist'] = thrustlist
    lists['isplist'] = isplist
    lists['rocketmachlist'] = rocketmachlist
    lists['alist'] = alist
    lists['timelist'] = timelist
    lists['burntimelist'] = burntimelist
    lists['xlist'] = xlist
    lists['rlist'] = rlist
    lists['machlist'] = machlist
    lists['preslist'] = preslist
    lists['templist'] = templist
    lists['chlist']  = chlist# (TC.rt/rlistflipped)**.5*.003 
    lists['twlist']  = twlist#(rlistflipped/TC.rt)*.001 
    lists['nlist']  = nlist#np.ones(len(xlist))*80
    lists['ewlist']  = ewlist#np.ones(len(xlist))*.005
    #HELIC'ITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    lists['helicitylist']  = helicitylist
    lists['alistflipped'] = alistflipped
    lists['vlistflipped'] = vlistflipped
    lists['hydraulicdiamlist'] = hydraulicdiamlist
    lists['salistflipped'] = salistflipped
    lists['dxlist'] = dxlist
    lists['cwlist'] = cwlist
    lists['Twglist'] = Twglist
    lists['hglist'] = hglist
    lists['qdotlist'] = qdotlist
    lists['Twclist'] = Twclist
    lists['hclist'] = hclist
    lists['Tclist'] = Tclist
    lists['coolantpressurelist'] = coolantpressurelist
    lists['fincoolingfactorlist'] = fincoolingfactorlist
    lists['rholist']=rholist
    lists['viscositylist']=viscositylist
    lists['Relist']=Relist
    lists['FOSlist']=FOSlist

    params['mi'] = args['drymass']
    params['L'] = L
    params['twg_max'] = np.max(Twglist)
    params['twc_max'] =  np.max(Twclist)
    return params, lists
