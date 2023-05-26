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

def main(args, configtitle, output=True):


    # FIRST DETERMINE INITIAL ESTIMATES FOR IDEAL PARAMS
    path=os.path.join( "Configs",configtitle)
    if output:
        os.makedirs(path,exist_ok=False)
    ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=output)
    if output:
        plt.savefig(os.path.join(path, "idealisp.png"))
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
    
    ### NEXT SOLVE FOR MIN IMPULSE TO REACH ALTITUDE GOAL
    miinit,lambdainit,totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank  = SA.mass_approx(params['pc'],params['dp'], 12, params['rho_fuel'],params['rho_ox'], params['thrust'], params['isp'], params['time'], params['rm'])
    #os.mkdir("Configs")
    #with open(os.path.join("..", "Configs","params.json"),"w") as file:
    #    json.dump({name: str(params[name]) for name in params}, file)
    #print(params)
    thrusttoweight_approx = args['TWR']

    #json.load()
    dt=.025

    params['thrust'] = totalmass*thrusttoweight_approx*9.81 #recompute with a resonable thrust
    params = FAC.SpreadsheetSolver(params)
    try: #if you inputted an impulseguess
        L, mi, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA_MassAprox(params, args['impulseguess'], thrusttoweight_approx, H=100000, dp = params['dp'],
                                                                            dt=dt, Af=None,
                                                                            ispcorrection=1)
        
        #impulselist = np.trim_zeros(impulselist)
        #hslist = np.trim_zeros(hslist)
        #iterations = np.arange(len(hslist))
        #plt.figure()
        #plt.plot(iterations.T,hslist,'g', label = "T_wg [K]")
        #plt.xlabel("Iteration")
        #plt.ylabel("Altitude Reached [m]")
        #plt.title("Altitude Convergance")
        #hslist.sort()
        #impulselist.sort()
        #plt.figure()
        #plt.plot(impulselist,hslist,'r*', label = "T_wg [K]")
        #plt.xlabel("Impulse [N*s]")
        #plt.ylabel("Altitude [m]")
        #plt.title("Recreation of Altitude vs Impulse from Iterations")

        #plt.show()
        
    except:
        L, mi, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA_MassAprox(params, 491801.7683699908, thrusttoweight_approx, H=100000, dp = params['dp'],
                                                                            dt=dt, Af=None,
                                                                            ispcorrection=1)
      
  # np.savetxt(os.path.join( path,"vlist.csv"), vlist, delimiter=",")
  # np.savetxt(os.path.join( path,"hlist.csv"), hlist, delimiter=",")                                                                              
    #plt.plot(impulselist,hslist)
    #plt.show()
    #Note that only thrust and time were changed by the rocket equation, so you need to reupdate params
    newargs = {
        'thrust': params['thrust'],  # Newtons
        'time': params['time'],  # s
        'pc': params['pc'],
        'pe': params['pe'],
        'rm' : params['rm'],
        'rc': params['rc'],
        'lstar': params['lstar'],
        'fuelname': params['fuelname'],
        'oxname': params['oxname'],
        'throat_radius_curvature': params['throat_radius_curvature'],
        'dp': params['dp']}
    params = FAC.SpreadsheetSolver(newargs)
    # now do it again with cr instead of rc to get all the specific values for finite combustors
    newargs = {
        'thrust': params['thrust'],  # Newtons
        'time': params['time'],  # s
        'pc': params['pc'],
        'pe': params['pe'],
        'rm' : params['rm'],
        'cr': params['cr'],
        'lstar': params['lstar'],
        'fuelname': params['fuelname'],
        'oxname': params['oxname'],
        'throat_radius_curvature': params['throat_radius_curvature'],
        'dp': params['dp']}
    params = FAC.SpreadsheetSolver(newargs)

    if output:
        with open(os.path.join(path, "mass_output.txt"), "w") as f:
            print(f"This is supposed to be a CSV! Eventually it should be passed params and get all that data too!", file=f)
             #overwriting current file
    mis, lambdas, totalmasses, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank = \
        SA.mass_approx(params['pc'], params['dp'], 12, params['rho_fuel'], params['rho_ox'],
                                                   params['thrust'], params['isp'], params['time'], params['rm'],
                   printoutput=output, outputdir=path)
    if output:
        with open(os.path.join(path, "mass_output.txt"), "a") as f:
            for param in list(params):
                if params[param] is None:
                    print(param + f", NONE", file=f)
                else:
                    #try:
                    #    print(param + f", " +'%.3f'%(params[param]), file=f)
                    #except:
                    print(param + f", {params[param]}", file=f)
            print(f"Deltav, {params['isp']*9.81*math.log(1/L)}", file=f)

        title = f"Fuel = {params['fuelname']}, " \
                f"Ox = {params['oxname']}, " \
                f"mi = {int(((1 / L) * params['mdot'] * params['time'] - params['mdot'] * params['time']))}," \
                f" mp = {int(params['mdot'] * params['time'])}, " \
                f"Mtotal = {int((1 / L) * params['mdot'] * params['time'])}"
        RE.ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time=params['time'], title=title, dt=dt)
        plt.savefig(os.path.join(path, "trajectory.png"))
    ### Now that we have the "optimal rocket", figure out flow and wall temps
    
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
    chlist  = (TC.rt/rlistflipped)**.5*.003 
    twlist  = (rlistflipped/TC.rt)*.001 
    nlist  = np.ones(len(xlist))*80
    ewlist  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylist  = (rlistflipped**1.5/TC.rt**1.5)*45*math.pi/180
    chanelToLandRatio = 2
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
    
    if output:
        #Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist = CS.steadyStateTemperatures(None,TC, params, salistflipped,n, coolingfactorlist,
        #                        heatingfactorlist, xlistflipped, vlistflipped ,293, params['pc']+params['dp'][0], twlistflipped, hydraulicdiamlist)

        title="ChamberTemps"
        fig, axs = plt.subplots(3,3)
        fig.suptitle(title)

        axs[0,1].plot(xlistflipped,hglist , 'g')  # row=0, column=0
        axs[1,1].plot(xlistflipped,hclist , 'r')# row=1, column=0
        axs[2,1].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0


        axs[0,1].set_title('hglist')
        axs[1,1].set_title('hclist')
        axs[2,1].set_title('Thrust Chamber Shape')

        axs[0,0].plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp")
        axs[0,0].plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp") # row=0, column=0
        axs[0,0].plot(xlistflipped,Tclist , 'b', label="Coolant Temp") #
        axs[1,0].plot(xlistflipped,Tclist , 'r')# row=1, column=0
        axs[2,0].plot(xlistflipped,hydraulicdiamlist , 'r')# row=1, column=0

        axs[0,0].set_title('Twg')
        axs[1,0].set_title('Tc')
        axs[2,0].set_title('hydraulicdiam')
        axs[0,0].legend()

        axs[0,2].plot(xlistflipped,Twglist*const.degKtoR-458.67 , 'g', label="Gas Side Wall Temp, F")
        #axs[0,2].plot(xlistflipped,Tcoatinglist*const.degKtoR-458.67 , 'k', label="Opposite side of coating Temp, F")
        axs[0,2].plot(xlistflipped,Twclist * const.degKtoR-458.67, 'r', label="CoolantSide Wall Temp, F") # row=0, column=0

        axs[0,2].plot(xlistflipped,Tclist * const.degKtoR-458.67, 'b', label="Coolant Temp, F") #
        axs[1,2].plot(xlistflipped,coolantpressurelist /const.psiToPa, 'k')
        axs[2,2].plot(xlistflipped,rholist, 'k') # row=0, column=0

        axs[0,2].set_title('Twg')
        axs[1,2].set_title('coolantpressure (psi)')
        axs[2,2].set_title('density of coolant')
        axs[0,2].legend()
        plt.savefig(os.path.join(path, "temperatures.png"))
        print(f"max twg = {np.max(Twglist)} in kelvin, {np.max(Twglist)*const.degKtoR} in Rankine (freedom)\n max Twc ="
              f" {np.max(Twclist)} in kelvin, {np.max(Twclist)*const.degKtoR} in Rankine (freedom)")
        # Hide x labels and tick labels for top plots and y ticks for right plots.

        title="Flow properties along thrust chamber"
        fig1, axs1 = plt.subplots(4,1)

        fig1.suptitle(title)

        axs1[0].plot(TC.xlist,machlist , 'g')  # row=0, column=0
        axs1[1].plot(TC.xlist,preslist , 'r')# row=1, column=0
        axs1[2].plot(TC.xlist,templist , 'b') # row=0, column=0
        axs1[3].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0


        axs1[0].set_title('Mach')
        axs1[1].set_title('Pressure')
        axs1[2].set_title('temperature')
        plt.savefig(os.path.join(path, "flowprops.png"))

        title=f"Chamber Wall Temperatures: Temp At Injector Face = {Twglist[-1]}"
        plt.figure()
        plt.plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp, K")
        plt.plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp, K") # row=0, column=0
        plt.plot(xlistflipped,Tclist , 'b', label="Coolant Temp, K") # 
        plt.xlabel("Axial Position [m From Injector Face]")
        plt.ylabel("Temperature [K]")
        plt.legend()
        plt.title(title)
        plt.savefig(os.path.join(path, "ChamberTemps_LachlanFormat.png"))

        axs[0,1].plot(xlistflipped,hglist , 'g')  # row=0, column=0
        axs[1,1].plot(xlistflipped,hclist , 'r')# row=1, column=0
        axs[2,1].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0


        axs[0,1].set_title('hglist')
        axs[1,1].set_title('hclist')
        axs[2,1].set_title('Thrust Chamber Shape')

        axs[0,0].plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp")
        axs[0,0].plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp") # row=0, column=0
        axs[0,0].plot(xlistflipped,Tclist , 'b', label="Coolant Temp") #
        axs[1,0].plot(xlistflipped,Tclist , 'r')# row=1, column=0
        axs[2,0].plot(xlistflipped,hydraulicdiamlist , 'r')# row=1, column=0

        axs[0,0].set_title('Twg')
        axs[1,0].set_title('Tc')
        axs[2,0].set_title('hydraulicdiam')
        axs[0,0].legend()

        axs[0,2].plot(xlistflipped,Twglist*const.degKtoR-458.67 , 'g', label="Gas Side Wall Temp, F")
        #axs[0,2].plot(xlistflipped,Tcoatinglist*const.degKtoR-458.67 , 'k', label="Opposite side of coating Temp, F")
        axs[0,2].plot(xlistflipped,Twclist * const.degKtoR-458.67, 'r', label="CoolantSide Wall Temp, F") # row=0, column=0

        axs[0,2].plot(xlistflipped,Tclist * const.degKtoR-458.67, 'b', label="Coolant Temp, F") #
        axs[1,2].plot(xlistflipped,coolantpressurelist /const.psiToPa, 'k')
        axs[2,2].plot(xlistflipped,rholist, 'k') # row=0, column=0

        axs[0,2].set_title('Twg')
        axs[1,2].set_title('coolantpressure (psi)')
        axs[2,2].set_title('density of coolant')
        axs[0,2].legend()
        plt.savefig(os.path.join(path, "temperatures.png"))
        print(f"max twg = {np.max(Twglist)} in kelvin, {np.max(Twglist)*const.degKtoR} in Rankine (freedom)\n max Twc ="
              f" {np.max(Twclist)} in kelvin, {np.max(Twclist)*const.degKtoR} in Rankine (freedom)")


        title=f"Chamber Wall Temperatures: Temp At Injector Face = {Twglist[-1]}"
        plt.figure()
        plt.plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp, K")
        plt.plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp, K") # row=0, column=0
        plt.plot(xlistflipped,Tclist , 'b', label="Coolant Temp, K") # 
        plt.xlabel("Axial Position [m From Injector Face]")
        plt.ylabel("Wall Temperature [K]")
        plt.title(title)
        plt.savefig(os.path.join(path, "ChamberTemps_LachlanFormat.png"))

        axs[0,1].plot(xlistflipped,hglist , 'g')  # row=0, column=0
        axs[1,1].plot(xlistflipped,hclist , 'r')# row=1, column=0
        axs[2,1].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0


        axs[0,1].set_title('hglist')
        axs[1,1].set_title('hclist')
        axs[2,1].set_title('Thrust Chamber Shape')

        axs[0,0].plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp")
        axs[0,0].plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp") # row=0, column=0
        axs[0,0].plot(xlistflipped,Tclist , 'b', label="Coolant Temp") #
        axs[1,0].plot(xlistflipped,Tclist , 'r')# row=1, column=0
        axs[2,0].plot(xlistflipped,hydraulicdiamlist , 'r')# row=1, column=0

        axs[0,0].set_title('Twg')
        axs[1,0].set_title('Tc')
        axs[2,0].set_title('hydraulicdiam')
        axs[0,0].legend()

        axs[0,2].plot(xlistflipped,Twglist*const.degKtoR-458.67 , 'g', label="Gas Side Wall Temp, F")
        #axs[0,2].plot(xlistflipped,Tcoatinglist*const.degKtoR-458.67 , 'k', label="Opposite side of coating Temp, F")
        axs[0,2].plot(xlistflipped,Twclist * const.degKtoR-458.67, 'r', label="CoolantSide Wall Temp, F") # row=0, column=0

        axs[0,2].plot(xlistflipped,Tclist * const.degKtoR-458.67, 'b', label="Coolant Temp, F") #
        axs[1,2].plot(xlistflipped,coolantpressurelist /const.psiToPa, 'k')
        axs[2,2].plot(xlistflipped,rholist, 'k') # row=0, column=0

        axs[0,2].set_title('Twg')
        axs[1,2].set_title('coolantpressure (psi)')
        axs[2,2].set_title('density of coolant')
        axs[0,2].legend()
        plt.savefig(os.path.join(path, "temperatures.png"))
        print(f"max twg = {np.max(Twglist)} in kelvin, {np.max(Twglist)*const.degKtoR} in Rankine (freedom)\n max Twc ="
              f" {np.max(Twclist)} in kelvin, {np.max(Twclist)*const.degKtoR} in Rankine (freedom)")

        plt.show()

    if ~output:
        params['mi'] = mis
        params['L'] = lambdas
        params['M'] = totalmasses
        params['wstruct'] = wstruct
        params['newheight'] = newheight
        params['heightox'] = heightox
        params['heightfuel'] = heightfuel
        params['vol_ox'] = vol_ox
        params['vol_fuel'] = vol_fuel
        params['P_tank'] = P_tank
        params['twg_max'] = np.max(Twglist)
        params['twc_max'] =  np.max(Twclist)
        return params

"""
args = {
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
       # 'phi':1,
        'cr': None,
        'TWR' : 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa,
        'impulseguess' : 800000,
        'rc' : .11} #623919}
configtitle = "NitroEthanol75 at 300 psi - 150 psi dp, 22 cm chamber diam"
main(args, configtitle, output=True)
"""

"""
mi = 10000
params = FAC.SpreadsheetSolver(args)
miinittemp = 0
pclist = np.linspace(450,1000,20)*const.psiToPa
milist = np.zeros(len(pclist))
isplist = np.zeros(len(pclist))
index=0
for pc in pclist:
    args = {
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': (pc-100*const.psiToPa)/1.2,
        'pe': 14.7 * const.psiToPa,
       # 'phi':1,
        'cr': 4,
        'TWR' : 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': [0,0],
        'impulseguess' : 800000} 
    ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=False)
    args['rm'] = mrideal
    params = FAC.SpreadsheetSolver(args)
    miinittemp = 0
    mi, L, totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank = SA.mass_approx(params['pc'], params['dp'], 12, params['rho_fuel'], params['rho_ox'],
                                        params['thrust'], params['isp'], params['time'], params['rm'])

    print(f"MI is {mi}, L is {L} and Thrust is {params['thrust']} and dp is {params['dp']}")
    milist[index] = mi
    isplist[index] = params['isp']
    index = index +1
plt.plot(pclist/const.psiToPa,milist)
plt.xlabel('Tank Pressure (psi)')
plt.ylabel('Dry Mass [kg]')
plt.title(f'Dry Mass vs Tank Pressure for NitroEthanol, 4000 lb thrust, 40s burn') #100 psi dp + 20\% injector

plt.figure()
plt.plot(pclist/const.psiToPa,isplist)
plt.xlabel('Tank Pressure (psi)')
plt.ylabel('isp [s]')
plt.title(f'isp vs Tank Pressure for NitroEthanol, 4000 lb thrust, 40s burn') #100 psi dp + 20\% injector
plt.show()
"""
"""configs = ethalox, ethanol nitrous, 
for pressure in np.arrange(100,700,25):
    for configs in configs
        main(pressure)
        
with open configs/massoutput.txt
    jon.load(massoutput.json)"""

"""
#RUNNING ROCKET EQUATION AND DETERMINING ALTIDUED FOR FIXED PARAMS

args = {
    'thrust': 1200 * const.lbToN,  # Newtons
    'time': 7.5,  # s
    # 'rho_ox' : 1141, #Kg/M^3
    # 'rho_fuel' : 842,
    'pc': 700 * const.psiToPa,
    'pe': 14.7 * const.psiToPa,
    'rm': 3.76/.95,
    # TYPICALLY WE SET 'phi'
    'cr': 4,
    'lstar': 1,
    'fuelname': 'Ethanol',
    'oxname': 'LOX',
    'throat_radius_curvature': .0254 / 2,
    'dp': 150 * const.psiToPa}

params = FAC.SpreadsheetSolver(args)
dp=250*const.psiToPa
miinit,lambdainit,totalmass = SA.mass_approx(params['pc'],dp, 12, params['rho_fuel'],params['rho_ox'], params['thrust'], params['isp'], params['time'], params['rm'])
#os.mkdir("Configs")
#with open(os.path.join("..", "Configs","params.json"),"w") as file:
#    json.dump({name: str(params[name]) for name in params}, file)
#print(params)
thrusttoweight_approx = 4

#json.load()
dt=.25

#params['thrust'] = totalmass*thrusttoweight_approx*9.81 #recompute with a resonable thrust
params = FAC.SpreadsheetSolver(params)
miinit, lambdainit, totalmass = SA.mass_approx(params['pc'], dp, 8, params['rho_fuel'], params['rho_ox'],
                                               params['thrust'], params['isp'], params['time'], params['rm'])
L, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA(params, mi = miinit ,thrust = params['thrust'], burntime = params['time'], L = lambdainit, H = None,
                                                                      dt=dt, Af=None,
                                                                      ispcorrection=.95)

L=lambdainit
title = f"Fuel = {params['fuelname']}, Thrust = {params['thrust']}, burntime = {params['time']}, " \
        f"mi = {((1 / L) * params['mdot'] * params['time'] - params['mdot'] * params['time'])}, mp = {params['mdot'] * params['time']}, Mtotal = {(1 / L) * params['mdot'] * params['time']}"
RE.ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time=params['time'], title=title, dt=dt)
"""

""" RUNNING FLOW SIMPLE AND GETTING THE PLOTS"""
"""
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
        'fuelname': 'Ethanol',
        'oxname': 'LOX',
        'throat_radius_curvature': .02}

params = FAC.SpreadsheetSolver(args)

xlist = np.linspace(0, params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], 1000)

rlist = RListGenerator.sharpRList(xlist, params['lc'], params['rc'], params['lc']+(params['rc']-params['rt'])/(2**.5),
                                  params['rt'], params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], params['re']) #xlist, xns, rc, xt, rt, xe, re
TC = ThrustChamber.ThrustChamber(rlist,xlist)
print(TC.rt, TC.xt, TC.xns)

machlist,preslist,templist = TC.flowSimple(params)

print(machlist[499]-machlist[498])

title="testing flowSimple for thrustchamber"
fig, axs = plt.subplots(4,1)
fig.suptitle(title)

axs[0].plot(TC.xlist,machlist , 'g')  # row=0, column=0
axs[1].plot(TC.xlist,preslist , 'r')# row=1, column=0
axs[2].plot(TC.xlist,templist , 'b') # row=0, column=0
axs[3].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0


axs[0].set_title('Mach')
axs[1].set_title('Pressure')
axs[2].set_title('temperature')
plt.show()

for ax in axs.flat:
    ax.set(xlabel='x')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
"""
"""
main()
"""
"""
#print(f"{[Twgi, hgi, hci, nf, Qdotguess, Qdothc, abs(Qdothc - Qdotguess)]},")
guesses = np.array([[1692.8250165869413, 1210.4187502907312, 12272.708286129538, 0.013781435608251378, 25413.022206704107, 219334.1805378995, 193921.1583311954],
[2275.5229222813546, 1116.4617684414663, 12272.708286129538, 0.013781435608251378, 14486.947591283237, 325355.25436037, 310868.30676908675],
[1332.697905694413, 1279.122651686839, 12272.708286129538, 0.013781435382462546, 33195.194572006505, 152844.70569451494, 119649.51112250844],
[1110.1271108925282, 1326.681677481956, 12272.708286129538, 0.013781431856309942, 38493.266503978004, 106827.23001273624, 68333.96350875823],
[972.570794801885, 1358.3407770893996, 12272.708286129538, 0.013781426005795017, 41983.37339849354, 72201.18575578938, 30217.812357295843],
[887.5563160906433, 1378.8625226610375, 12272.708286129538, 0.013781422125123902, 44230.958609416615, 49584.528528244904, 5353.569918828289],
[835.0144787112416, 1391.934340905879, 11440.141117134803, 0.013781422225300682, 45656.802748091875, 32998.14535571232, 12658.657392379559],
[883.4461991982033, 1379.8741446642316, 12272.708286129538, 0.013781421937719902, 44341.46335403265, 48473.99123975513, 4132.527885722484],
[869.0054142575317, 1383.4430413664822, 12272.708286129538, 0.013781421286061322, 44731.09798728304, 44569.67201492076, 161.42597236228175],
[856.0220321882558, 1386.6712917950035, 12272.708286129538, 0.013781420741277802, 45083.25587111174, 41106.44171831652, 3976.8141527952175],
[869.5972507676269, 1383.2963265637009, 12272.708286129538, 0.013781421310417069, 44715.08695726494, 44726.433689269455, 11.346732004516525],
[874.8870783596626, 1381.9866965160932, 12272.708286129538, 0.013781421536487887, 44572.14170594803, 46139.13389666865, 1566.9921907206226],
[870.6630102427988, 1383.0322243449277, 12272.708286129538, 0.013781421354173408, 44686.2639860794, 45008.48687602289, 322.2228899434922],
[869.6867670127857, 1383.2741390552292, 12272.708286129538, 0.013781421314097365, 44712.66557891686, 44750.13592069225, 37.47034177539172],
[869.4597252436268, 1383.3304153715214, 12272.708286129538, 0.013781421304761142, 44718.80712914407, 44690.01528649571, 28.791842648359307],
[869.5628820315197, 1383.3048454368134, 12272.708286129538, 0.013781421309003814, 44716.01663979802, 44717.33291522701, 1.3162754289878649],
[869.5441143385166, 1383.3094973819782, 12272.708286129538, 0.013781421308232022, 44716.52431555377, 44712.36313238565, 4.161183168114803],
[869.5595358004846, 1383.3056748637518, 12272.708286129538, 0.013781421308866209, 44716.1071568023, 44716.44682278205, 0.3396659797435859],
[869.5561340248978, 1383.306518059722, 12272.708286129538, 0.013781421308726317, 44716.1991764299, 44715.54601887051, 0.6531575593908201]])
iterations = np.arange(guesses.shape[0])

plt.figure()
plt.plot(iterations.T,guesses[:,0],'g', label = "T_wg [K]")
plt.plot(iterations.T,guesses[:,1],'r', label = "hg [w/(m^2 K)]")
plt.xlabel("Iteration")
plt.legend()
plt.title("Temperature and Heat Transfer Convergance")

plt.figure()
plt.plot(iterations.T,guesses[:,-1],'g', label = "T_wg [K]")
plt.xlabel("Iteration")
plt.ylabel("Error [W/m]")
plt.yscale('log')
plt.title("Heat Flow Error Convergance : abs(Qdothc - Qdotguess) ")

plt.show()
"""