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
def main(args, configtitle, output=True):


    # FIRST DETERMINE INITIAL ESTIMATES FOR IDEAL PARAMS
    if output:
        path=os.path.join( "Configs",configtitle)
        os.makedirs(path,exist_ok=False)
    ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen = DOMR.optimalMr(args, plot=output)
    plt.savefig(os.path.join(path, "idealisp.png"))
    print(f"isp max = {ispmaxavg}, ideal mr is {mrideal}")
    args['rm']=mrideal
    params = FAC.SpreadsheetSolver(args)

    ### NEXT SOLVE FOR MIN IMPULSE TO REACH ALTITUDE GOAL
    miinit,lambdainit,totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank  = SA.mass_approx(params['pc'],params['dp'], 12, params['rho_fuel'],params['rho_ox'], params['thrust'], params['isp'], params['time'], params['rm'])
    #os.mkdir("Configs")
    #with open(os.path.join("..", "Configs","params.json"),"w") as file:
    #    json.dump({name: str(params[name]) for name in params}, file)
    #print(params)
    thrusttoweight_approx = args['TWR']

    #json.load()
    dt=.1

    params['thrust'] = totalmass*thrusttoweight_approx*9.81 #recompute with a resonable thrust
    params = FAC.SpreadsheetSolver(params)
    try: #if you inputted an impulseguess
        L, mi, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA_MassAprox(params, args['impulseguess'], thrusttoweight_approx, H=100000, dp = params['dp'],
                                                                            dt=dt, Af=None,
                                                                            ispcorrection=.95)
    except:
        L, mi, hlist, vlist, thrustlist, isplist, machlist = RE.rocketEquationCEA_MassAprox(params, 491801.7683699908, thrusttoweight_approx, H=100000, dp = params['dp'],
                                                                            dt=dt, Af=None,
                                                                            ispcorrection=.95)
                                                                                    
    #plt.plot(impulselist,hslist)
    #plt.show()
    #Note that only thrust and time were changed by the rocket equation, so you need to reupdate params
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
    xlist = np.linspace(0, params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], 1000)

    rlist = RListGenerator.sharpRList(xlist, params['lc'], params['rc'], params['lc']+(params['rc']-params['rt'])/(2**.5),
                                      params['rt'], params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], params['re']) #xlist, xns, rc, xt, rt, xe, re
    TC = ThrustChamber.ThrustChamber(rlist,xlist)
    print(TC.rt, TC.xt, TC.xns)

    machlist,preslist,templist = TC.flowSimple(params)







    chanelthickness = (rlist==TC.rc)*.05*.0254+(rlist!=TC.rc)*.05*.0254
    tw=.025*.0254
    alist=math.pi*((rlist+tw+chanelthickness)**2-(rlist+tw)**2)
    alistflipped, n, coolingfactorlist, heatingfactorlist, xlistflipped, vlistflipped, twlistflipped, hydraulicdiamlist, salistflipped = CS.coaxialShellSetup(TC,
                                            params,rlist, tw=tw, chanelthickness = chanelthickness, helicity=1, dt=.005,vlist = None, alist=alist)

    Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist = CS.steadyStateTemperatures(None,TC, params, salistflipped,n, coolingfactorlist,
                                heatingfactorlist, xlistflipped, vlistflipped ,293, params['pc']+params['dp'], twlistflipped, hydraulicdiamlist)

    if output:
        title="testing flowSimple for thrustchamber"
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
args = {
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
       # 'phi':1,
        'cr': 4,
        'TWR' : 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa,
        'impulseguess' : 623919}
configtitle = "NitroEthanol75 at 300 psi"
main(args, configtitle, output=True)

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

