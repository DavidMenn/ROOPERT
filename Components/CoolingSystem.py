"""Here will lie the wall temps. Early on we're just going to do throat temp just to get sizing down"""
import sys
sys.path.insert(1,"./")
import json
""" Cooldown analysis conops
Use thrust chamber to get flow props and store geometry
pass chamber props to setup function for whatever manufacturing technique
acquire cooling passage properties
send them to the trenches of steady state cooling
return heats
send them to the structural analysis
return margin of safety
iterate until we reach the ideal thickness
"""
import Toolbox.Constant as const
import scipy
import re as regex
from rocketprops.rocket_prop import get_prop
import numpy as np
import math
from scipy import interpolate
import Toolbox.IsentropicEquations as Ise
import matplotlib.pyplot as plt


def coaxialShellSetup(thrustchamber, params, rlist, tw, chanelthickness, helicity, dt, vlist=None, alist=None):
    if vlist is None:
        vlist = params['mdot_fuel'] / params['rho_fuel'] / alist
    else:
        alist = params['mdot_fuel'] / params['rho_fuel'] / vlist
    n = 1

    twlist = tw * np.ones((1, rlist.size))
    vInterpolator = interpolate.interp1d(thrustchamber.xlist,
                                         vlist, kind='linear')
    rInterpolator = interpolate.interp1d(thrustchamber.xlist,
                                         rlist, kind='linear')
    twInterpolator = interpolate.interp1d(thrustchamber.xlist,
                                          twlist, kind='linear')
    aInterpolator = interpolate.interp1d(thrustchamber.xlist,
                                         alist, kind='linear')
    chanelthicknessInterpolator = interpolate.interp1d(thrustchamber.xlist,
                                         chanelthickness, kind='linear')
    x = thrustchamber.xlist[-1]
    ind = 0
    #xlist = np.zeros(int(thrustchamber.xlist[-1] / dt / np.min(vlist)))
    #dxlist = np.zeros(int(thrustchamber.xlist[-1] / dt / np.min(vlist)))
    xlist = np.zeros(int(thrustchamber.xlist[-1] / dt )+1)
    dxlist = np.zeros(int(thrustchamber.xlist[-1] / dt )+1)

    xlist[0] = x
    dx = dt
    while x > 0 + dx*1.01:

        dydx = (rInterpolator(x) - rInterpolator(x - dx)) / dx
        dhypotenuse = vInterpolator(x) * dt
        dx = dt#(dhypotenuse ** 2 / ((dydx ** -2) + 1) / (dydx ** 2)) ** .5
        if dydx == 0:
            dx = dt#dhypotenuse

        dxlist[ind] = dx
        x = xlist[ind] - dx
        xlist[ind + 1] = x
        ind = ind + 1

    while xlist[-1] == 0 and xlist[-2] == 0:
        xlist = xlist[0:-1]
        dxlist = dxlist[0:-1]  # trims the extra zeros
    dxlist[-1] = dxlist[-2]
    alistflipped = aInterpolator(xlist)
    vlistflipped = vInterpolator(xlist)
    xlistflipped = xlist
    twlistflipped = twInterpolator(xlist)[0]
    salistflipped = math.pi * 2 * (rInterpolator(xlist) + twlistflipped) * dxlist
    hydraulicdiamlist = 2 * (chanelthicknessInterpolator(xlist)) # Hydraulic diam is Douter-Dinner for an anulus by wetted perimiter
    coolingfactorlist = np.ones(xlist.size)
    heatingfactorlist = np.ones(xlist.size)   # .6 is from cfd last year, i think its bs but whatever

    return alistflipped, n, coolingfactorlist, heatingfactorlist, xlistflipped, vlistflipped, twlistflipped, hydraulicdiamlist, salistflipped

def finApproximation(geometry):

    return coolingfactorlist
# def brazedTubingSetup():
#    return alist, n, coolingfactorlist, heatingfactorlist, xlist, dxlist, tlist

# def milledChanelSetup():
#    return alist, n, coolingfactorlist, heatingfactorlist, xlist, dxlist

def steadyStateTemperatures(wallmaterial, thrustchamber, params, salist, numchanels, coolingfactorlist,
                            heatingfactorlist, xlist, vlist, initialcoolanttemp, initialcoolantpressure, twlist,
                            hydraulicdiamlist):
    if not (get_prop(params['fuelname']) is None):
        pObj = get_prop(params['fuelname'])
        ethpercent = None
        pObjWater = None
    else:
        try:
            print('Fuel does not exist, assuming its a water ethanol blend')
            ethpercent = int(regex.search(r'\d+', params['fuelname']).group()) / 100
            pObj = get_prop("ethanol")
            pObjWater = get_prop("water")
        except:
            print("rocketProps busted lol")

    try:
        wallmaterialprops = json.load(wallmaterial)
    except:
        print("either you didn't call this with a mat name or its not set up yet")
        kwInterpolator = interpolate.interp1d(np.hstack((25,np.arange(100, 1001, 100) + 273.5, np.array([10000]))),
                                              np.array([12.9, 13.9, 16.1, 18.2, 20, 22.1, 24.7, 31.4, 30.1, 31.5, 36.7,
                                                        36.7]), kind='linear')
        matprops = {'kw': kwInterpolator}  # copper = 398, mild steel = 51??

    staticnozzleparameters = {
        'throatdiameter': thrustchamber.rt * 2,
        'prandtlns': params['prns'],
        'viscosityns': params['viscosityns'],
        'cpns': params['cpns'],
        'pcns': params['pc'],
        'cstar': params['cstar'],
        'throatRadiusCurvature': params["throat_radius_curvature"],
        'at': thrustchamber.at,
        'tcns': params['temp_c'],
        'kw': matprops['kw'],
        'pObj': pObj,
        'pObjWater': pObjWater,
        'ethpercent': ethpercent,
        'roughness' : 32*10**-6
    }

    Trlist = np.zeros(xlist.size)
    for ind in np.arange(0, xlist.size):
        x = xlist[ind]
        Trlist[ind] = recoveryTemp(thrustchamber.tempInterpolator(xlist[ind]), params['gamma'],
                                   thrustchamber.machInterpolator(xlist[ind]), Pr=params['pr_throat'])

    Tc = initialcoolanttemp
    coolantpressure = initialcoolantpressure
    coolantpressurelist = np.zeros(xlist.size)
    Twglist = np.zeros(xlist.size)
    hglist = np.zeros(xlist.size)
    qdotlist = np.zeros(xlist.size)
    Twclist = np.zeros(xlist.size)
    hclist = np.zeros(xlist.size)
    Tclist = np.zeros(xlist.size)
    rholist = np.zeros(xlist.size)
    viscositylist = np.zeros(xlist.size)
    Relist = np.zeros(xlist.size)
    """
    # THIS IS JUST TO output CF plot
    cfindex=0
    velocitieslist = np.arange(.5,100,.1)
    cflist= np.zeros(velocitieslist.size)
    reynoldslist = np.zeros(velocitieslist.size)
    for velocity in velocitieslist:
        cflist[cfindex]=cf(Tc, coolantpressure, velocity, hydraulicdiamlist[0],
           staticnozzleparameters, roughness=32*10**-6)
        reynoldslist[cfindex]=reynolds(Tc,coolantpressure,velocity,hydraulicdiamlist[0],pObj, pObjWater=None, ethpercent=None)
        cfindex=cfindex+1
    title = f"cf for relative roughenss ={32*10**-6/hydraulicdiamlist[0]}"
    plt.figure()

    plt.plot(reynoldslist,cflist*4, 'g')  # row=0, column=0
    plt.title(title)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True)
    plt.show()
    """
    for ind in np.arange(0, xlist.size):
        Tclist[ind] = Tc  # set cooland at current station to coolant temp
        if ind == 0:
            coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                       staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                      xlist[ind] - xlist[ind + 1]) / \
                              hydraulicdiamlist[ind] * (
                                          .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)
        else:
            coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                       staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                      xlist[ind - 1] - xlist[ind]) / \
                              hydraulicdiamlist[ind] * (
                                      .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)


        coolantpressurelist[ind] = coolantpressure
        x = xlist[ind]
        print(x)
        Tri = Trlist[ind]
        gassidearea = thrustchamber.areaInterpolator(x)

        Twglist[ind] = qdotdiffMinimizer(staticnozzleparameters, gassidearea, thrustchamber.machInterpolator(x),
                                         params['gamma'],
                                         Tri, Tc, twlist[ind], coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                         coolingfactorlist[ind],heatingfactorlist[ind])
        hglist[ind] = heatingfactorlist[ind]*Bartz(staticnozzleparameters['throatdiameter'], staticnozzleparameters['viscosityns'],
                            staticnozzleparameters['prandtlns'], staticnozzleparameters['cpns'],
                            staticnozzleparameters['pcns'], staticnozzleparameters['cstar'],
                            staticnozzleparameters['throatRadiusCurvature'], staticnozzleparameters['at'],
                            gassidearea, Twglist[ind], staticnozzleparameters['tcns'],
                            thrustchamber.machInterpolator(x), params['gamma'])
        qdotlist[ind] = hglist[ind] * (Trlist[ind] - Twglist[ind])
        # This works since its the same method used in qdotdiff function, its just a one iteration approx
        Twclist[ind] = Twglist[ind] - qdotlist[ind] / (matprops['kw']((Twglist[ind])) / twlist[ind])
        Twclist[ind] = Twglist[ind] - qdotlist[ind] / (matprops['kw']((Twglist[ind]+Twclist[ind])/2) / twlist[ind])
        hclist[ind] = hc(Twclist[ind], Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                         staticnozzleparameters['pObj'], staticnozzleparameters['pObjWater'],
                         staticnozzleparameters['ethpercent'], coolingfactorlist[ind])
        Tc = Tc + qdotlist[ind] * salist[ind] / (params['mdot_fuel'] * heatCapacity(Tc, staticnozzleparameters))
        rholist[ind] = rho(Tc, coolantpressure, staticnozzleparameters)
        viscositylist[ind]= viscosity(Tc,coolantpressure,  staticnozzleparameters['pObj'],staticnozzleparameters['pObjWater'],staticnozzleparameters['ethpercent'])
        Relist[ind] = reynolds(Tc,coolantpressure,vlist[ind],hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                                       staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])

    return Twglist, hglist, qdotlist, Twclist, hclist, Tclist, \
           coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist

def ablative():
    return "set this up"


def transientTemperature(wall, steadystatestuff):
    "solve steady state with twg set to room temp at first with some time step"


# ns stands for nozzle start, use finite area combustor CEA to find it (or just use pc and tcomb)
# This equation neglects g because you dont need it for si units (i think)
def Bartz(throatdiameter, viscosityns, prandtlns, cpns, pcns, cstar, throatRadiusCurvature, at, a, twg, tcns, mach,
          gamma):
    sigma = ((.5 * twg / tcns * (1 + (mach ** 2) * (gamma - 1) / 2) + .5) ** (-.68)) * (
            1 + (mach ** 2) * (gamma - 1) / 2) ** (-.12)
    return (.026 / (throatdiameter ** .2) * (cpns * viscosityns ** .2) / (prandtlns ** .6) * (pcns / cstar) ** (.8) * (
            throatdiameter / throatRadiusCurvature) ** .1) * (at / a) ** .9 * sigma


def recoveryTemp(temp, gam, mach, Pr=None):
    Taw = Ise.totalT(temp, gam, mach) / (1 + (((gam - 1) / 2) * (mach ** 2)))
    if Pr is None:
        Pr = .6865
    r = Pr ** (1 / 3)
    return Taw * (1 + (r * (mach ** 2) * ((gam - 1) / 2)))


def qdotdiffMinimizer(staticnozzleparams, a, mach, gamma,
                      Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam, coolingfactor, heatingfactor):
    return scipy.optimize.minimize_scalar(qdotdiff, bounds=(750, Tri), tol=10, method='bounded',
                                          args=(staticnozzleparams, a, mach, gamma,
                                                Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam,
                                                coolingfactor, heatingfactor))['x']  # going to find root of qdotdiff for twg
def qdotdiffMinimizerJayLeno(staticnozzleparams, a, mach, gamma,
                      Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam, coolingfactor, heatingfactor):
    return scipy.optimize.minimize_scalar(qdotdiffjayleno, bounds=(750, Tri), method='bounded', tol=10, args=(staticnozzleparams, a, mach, gamma,
                                                Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam,
                                                coolingfactor, heatingfactor))['x']

def qdotdiff(Twgi, staticnozzleparams, a, mach, gamma,
             Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam, coolingfactor, heatingfactor):
    hgi = heatingfactor*Bartz(staticnozzleparams['throatdiameter'], staticnozzleparams['viscosityns'],
                staticnozzleparams['prandtlns'], staticnozzleparams['cpns'],
                staticnozzleparams['pcns'], staticnozzleparams['cstar'], staticnozzleparams['throatRadiusCurvature'],
                staticnozzleparams['at'],
                a, Twgi, staticnozzleparams['tcns'], mach, gamma)
    qdotguess = hgi * (Tri - Twgi)
    Twci = Twgi - ((qdotguess * tw) / staticnozzleparams['kw'](Twgi))# get the initial guess
    Twci = Twgi - ((qdotguess * tw) / staticnozzleparams['kw']((Twgi+Twci)/2)) # now use avg wall temp
    hci = hc(Twci, Tc, coolantpressure, coolantvelocity, hydraulicdiam, staticnozzleparams['pObj'],
             staticnozzleparams['pObjWater'], staticnozzleparams['ethpercent'], coolingfactor)
    qdothc = hci * (Twci - Tc)
    return abs(qdothc - qdotguess)

#this is adding a layer of cermaic coating (similar to the jay leno steam fireb ox thing idk its in a video somewhere)
def qdotdiffjayleno(Twgi, staticnozzleparams, a, mach, gamma,
             Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam, coolingfactor, heatingfactor):
    hgi = heatingfactor*Bartz(staticnozzleparams['throatdiameter'], staticnozzleparams['viscosityns'],
                staticnozzleparams['prandtlns'], staticnozzleparams['cpns'],
                staticnozzleparams['pcns'], staticnozzleparams['cstar'], staticnozzleparams['throatRadiusCurvature'],
                staticnozzleparams['at'],
                a, Twgi, staticnozzleparams['tcns'], mach, gamma)
    hgi = 1/(1/hgi+1/staticnozzleparams['hcoating'])
    qdotguess = hgi * (Tri - Twgi)
    Twci = Twgi - ((qdotguess * tw) / staticnozzleparams['kw'](Twgi))  # get the initial guess
    Twci = Twgi - ((qdotguess * tw) / staticnozzleparams['kw']((Twgi + Twci) / 2))  # now use avg wall temp
    hci = hc(Twci, Tc, coolantpressure, coolantvelocity, hydraulicdiam, staticnozzleparams['pObj'],
             staticnozzleparams['pObjWater'], staticnozzleparams['ethpercent'], coolingfactor)
    qdothc = hci * (Twci - Tc)
    return abs(qdothc - qdotguess)
def hc(tempwall, temp, pres, fluidvelocity, hydraulicdiam, pObj, pObjWater, ethpercent, coolingfactor):
    Re = reynolds(temp, pres, fluidvelocity, hydraulicdiam, pObj, pObjWater, ethpercent)
    Pr = prandtl(temp, pres, pObj, pObjWater, ethpercent)
    if ethpercent is None:
        thermalcond = pObj.CondAtTdegR(temp * const.degKtoR) * const.BTUperHrFtRtoWperMK
        viscosity = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                              pres / const.psiToPa)
        viscosityw = .1 * pObj.Visc_compressed(tempwall * const.degKtoR,
                                               pres / const.psiToPa)
        if math.isnan(viscosityw):
            viscosityw=viscosity #make this a non factor since it keeps breaking
    else:
        thermalcondeth = pObj.CondAtTdegR(temp * const.degKtoR) * const.BTUperHrFtRtoWperMK
        thermalcondwater = pObjWater.CondAtTdegR(temp * const.degKtoR) * const.BTUperHrFtRtoWperMK
        thermalcond = (ethpercent * thermalcondeth + (1 - ethpercent) * thermalcondwater)
        viscosityeth = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                                 pres / const.psiToPa)
        viscositywater = .1 * pObjWater.Visc_compressed(temp * const.degKtoR,
                                                        pres / const.psiToPa)
        viscosityethw = .1 * pObj.Visc_compressed(tempwall * const.degKtoR,
                                                  pres / const.psiToPa)
        viscositywaterw = .1 * pObjWater.Visc_compressed(tempwall * const.degKtoR,
                                                         pres / const.psiToPa)
        viscosity = ethpercent * viscosityeth + (1 - ethpercent) * viscositywater
        viscosityw = ethpercent * viscosityethw + (1 - ethpercent) * viscositywaterw
        if math.isnan(viscosityw):
            viscosityw=viscosity #make this a non factor since it keeps breaking
    # FIGURE OUT THESE COEFFICIENTS PETE, PAGE 198 IN HEISTER eider - tate equation
    a = .027  # this is the weirtd one, its just a multiplier. Maybe include thisi n the cooling factor list? it depends on coolabnt, heister page 198
    m = .8
    n = .5 #.4
    b = .114
    #if pObj.PvapAtTdegR(temp * const.degKtoR) * const.psiToPa > pres:  # THIS MEANS YOUR ETHNAOL IS BOILING
    #    coolingfactor = coolingfactor / 10  # this is made up but hopefully will make it obvious your shit is boiling
    return coolingfactor * thermalcond / hydraulicdiam * a * Re ** m * Pr ** n * (viscosity / viscosityw) ** b


def reynolds(temp, pres, fluidvelocity, hydraulicdiameter, pObj, pObjWater, ethpercent):
    if ethpercent is None:
        rho = 1000 * pObj.SG_compressed(temp * const.degKtoR,
                                        pres / const.psiToPa)
        viscosity = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                              pres / const.psiToPa)
    else:
        rhoeth = 1000 * pObj.SG_compressed(temp * const.degKtoR,
                                           pres / const.psiToPa)
        viscosityeth = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                                 pres / const.psiToPa)
        rhowater = 1000 * pObjWater.SG_compressed(temp * const.degKtoR,
                                                  pres / const.psiToPa)
        viscositywater = .1 * pObjWater.Visc_compressed(temp * const.degKtoR,
                                                        pres / const.psiToPa)
        viscosity = ethpercent * viscosityeth + (1 - ethpercent) * viscositywater
        rho = ethpercent * rhoeth + (1 - ethpercent) * rhowater
    return rho * fluidvelocity * hydraulicdiameter / viscosity

# W/mK, kg/m3, J/kgK
def prandtl(temp, pres, pObj, pObjWater, ethpercent):
    if ethpercent is None:
        rho = 1000 * pObj.SG_compressed(temp * const.degKtoR,
                                        pres / const.psiToPa)
        thermalcond = pObj.CondAtTdegR(temp * const.degKtoR) * const.BTUperHrFtRtoWperMK
        cp = pObj.CpAtTdegR(temp * const.degKtoR) * const.BTUperLbmRtoJperKgK
        thermaldiffusivity = thermalcond / (rho * cp)  # this the definition
        viscosity = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                              pres / const.psiToPa)  # the .1 is because its giving it in units of poise not SI
    else:
        rhoeth = 1000 * pObj.SG_compressed(temp * const.degKtoR,
                                           pres / const.psiToPa)
        viscosityeth = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                                 pres / const.psiToPa)
        thermalcondeth = pObj.CondAtTdegR(temp * const.degKtoR) * const.BTUperHrFtRtoWperMK
        cpeth = pObj.CpAtTdegR(temp * const.degKtoR) * const.BTUperLbmRtoJperKgK
        rhowater = 1000 * pObjWater.SG_compressed(temp * const.degKtoR,
                                                  pres / const.psiToPa)
        viscositywater = .1 * pObjWater.Visc_compressed(temp * const.degKtoR,
                                                        pres / const.psiToPa)
        thermalcondwater = pObjWater.CondAtTdegR(temp * const.degKtoR) * const.BTUperHrFtRtoWperMK
        cpwater = pObjWater.CpAtTdegR(temp * const.degKtoR) * const.BTUperLbmRtoJperKgK
        thermaldiffusivity = (ethpercent * thermalcondeth + (1 - ethpercent) * thermalcondwater) / (
                    (ethpercent * rhoeth + (1 - ethpercent) * rhowater) * (
                        ethpercent * cpeth + (1 - ethpercent) * cpwater))  # this the definition
        viscosity = ethpercent * viscosityeth + (1 - ethpercent) * viscositywater
        rho = ethpercent * rhoeth + (1 - ethpercent) * rhowater
    return viscosity / thermaldiffusivity / rho  # viscosity over rho to go from dyunamic vsicosity to kinematic viscosity (momentum diffusivity)
    #return viscosity*((ethpercent * cpeth + (1 - ethpercent) * cpwater))/ (ethpercent * thermalcondeth + (1 - ethpercent) * thermalcondwater)

def heatCapacity(temp, staticnozzleparameters):
    if staticnozzleparameters['ethpercent'] is None:
        return staticnozzleparameters['pObj'].CpAtTdegR(temp * const.degKtoR) * const.BTUperLbmRtoJperKgK
    else:
        cpeth = staticnozzleparameters['pObj'].CpAtTdegR(temp * const.degKtoR) * const.BTUperLbmRtoJperKgK
        cpwater = staticnozzleparameters['pObjWater'].CpAtTdegR(temp * const.degKtoR) * const.BTUperLbmRtoJperKgK
        return (staticnozzleparameters['ethpercent'] * cpeth + (1 - staticnozzleparameters['ethpercent']) * cpwater)


def rho(temp, pres, staticnozzleparameters):
    if staticnozzleparameters['ethpercent'] is None:
        rho = 1000 * staticnozzleparameters['pObj'].SG_compressed(temp * const.degKtoR,
                                                                  pres / const.psiToPa)
        return rho
    else:
        rhoeth = 1000 * staticnozzleparameters['pObj'].SG_compressed(temp * const.degKtoR,
                                                                     pres / const.psiToPa)

        rhowater = 1000 * staticnozzleparameters['pObjWater'].SG_compressed(temp * const.degKtoR,
                                                                            pres / const.psiToPa)

        return staticnozzleparameters['ethpercent'] * rhoeth + (1 - staticnozzleparameters['ethpercent']) * rhowater


def cf(temp, pres, fluidvelocity, hydraulicdiameter, staticnozzleparameters, roughness):
    Re = reynolds(temp, pres, fluidvelocity, hydraulicdiameter, staticnozzleparameters['pObj'],
                  staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])
    if Re < 2100:
        return 16 / Re
    else:
        return scipy.optimize.minimize_scalar(turbulentCfImplicit, args=(roughness, hydraulicdiameter, Re))['x']
    # FIX THIS WHY IS IT GIVING WAY TOO HIGH PRESSURE DROP

def turbulentCfImplicit(cf, roughness, diameter,
                        Re):  # from wikipedia of moody chart lol, we are using fanning friction factor
    try:
        return abs(
            -2 * math.log(roughness / diameter / 3.7 + 2.51 / Re / math.sqrt(cf * 4), 10) - 1 / math.sqrt(cf * 4))
    except:
        return abs(
            -2 * math.log(roughness / diameter / 3.7 + 2.51 / Re / math.sqrt(cf * 4), 10) - 1 / math.sqrt(
                (cf + .00001) * 4))
def viscosity(    temp, pres,  pObj, pObjWater, ethpercent):
    if ethpercent is None:
        return .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                              pres / const.psiToPa)
    else:

        viscosityeth = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                                 pres / const.psiToPa)
        viscositywater = .1 * pObjWater.Visc_compressed(temp * const.degKtoR,
                                                        pres / const.psiToPa)

        return  ethpercent * viscosityeth + (1 - ethpercent) * viscositywater

# DELETE THIS SOON AND ADD THE COATING ANOTHER WAY
def steadyStateTemperaturesJayLeno(wallmaterial, thrustchamber, params, salist, numchanels, coolingfactorlist,
                            heatingfactorlist, xlist, vlist, initialcoolanttemp, initialcoolantpressure, twlist,
                            hydraulicdiamlist):
    if not (get_prop(params['fuelname']) is None):
        pObj = get_prop(params['fuelname'])
        ethpercent = None
        pObjWater = None
    else:
        try:
            print('Fuel does not exist, assuming its a water ethanol blend')
            ethpercent = int(regex.search(r'\d+', params['fuelname']).group()) / 100
            pObj = get_prop("ethanol")
            pObjWater = get_prop("water")
        except:
            print("rocketProps busted lol")

    try:
        wallmaterialprops = json.load(wallmaterial)
    except:
        #IT WOULD BE COOL TO HAVE LOADABLE MATERIALS!
        # print("either you didn't call this with a mat name or its not set up yet")
        kwInterpolator = interpolate.interp1d(np.hstack(np.arange(100,1001,100)+273.5,10000),
                                             np.array([12.9,13.9,16.1,18.2,20,22.1,24.7,31.4,30.1,31.5,36.7,36.7]), kind='linear')
        matprops = {'kw': kwInterpolator}  # copper = 398, mild steel = 51??

    staticnozzleparameters = {
        'throatdiameter': thrustchamber.rt * 2,
        'prandtlns': params['prns'],
        'viscosityns': params['viscosityns'],
        'cpns': params['cpns'],
        'pcns': params['pc'],
        'cstar': params['cstar'],
        'throatRadiusCurvature': params["throat_radius_curvature"],
        'at': thrustchamber.at,
        'tcns': params['temp_c'],
        'kw': matprops['kw'],
        'pObj': pObj,
        'pObjWater': pObjWater,
        'ethpercent': ethpercent,
        'coatingthickness' : .00076,
        'kcoating' : 10.1,
        'hcoating': 14470,
        'roughness' : 32*10**-6
    }

    Trlist = np.zeros(xlist.size)
    for ind in np.arange(0, xlist.size):
        x = xlist[ind]
        Trlist[ind] = recoveryTemp(thrustchamber.tempInterpolator(xlist[ind]), params['gamma'],
                                   thrustchamber.machInterpolator(xlist[ind]), Pr=params['pr_throat'])

    Tc = initialcoolanttemp
    coolantpressure = initialcoolantpressure
    coolantpressurelist = np.zeros(xlist.size)
    Twglist = np.zeros(xlist.size)
    hglist = np.zeros(xlist.size)
    qdotlist = np.zeros(xlist.size)
    Twclist = np.zeros(xlist.size)
    hclist = np.zeros(xlist.size)
    Tclist = np.zeros(xlist.size)
    rholist = np.zeros(xlist.size)
    viscositylist = np.zeros(xlist.size)
    Relist = np.zeros(xlist.size)

    for ind in np.arange(0, xlist.size):
        Tclist[ind] = Tc  # set cooland at current station to coolant temp
        if ind == 0:
            coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                       staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                      xlist[ind] - xlist[ind+1]) / \
                              hydraulicdiamlist[ind] * (.5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)
        else:
            coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                   staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                  xlist[ind - 1] - xlist[ind]) / \
                          hydraulicdiamlist[ind] * (
                                  .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)
        coolantpressurelist[ind] = coolantpressure
        x = xlist[ind]
        print(x)
        Tri = Trlist[ind]
        gassidearea = thrustchamber.areaInterpolator(x)

        Twglist[ind] = qdotdiffMinimizerJayLeno(staticnozzleparameters, gassidearea, thrustchamber.machInterpolator(x),
                                         params['gamma'],
                                         Tri, Tc, twlist[ind], coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                         coolingfactorlist[ind],heatingfactorlist[ind])
        hglist[ind] = heatingfactorlist[ind]*Bartz(staticnozzleparameters['throatdiameter'], staticnozzleparameters['viscosityns'],
                            staticnozzleparameters['prandtlns'], staticnozzleparameters['cpns'],
                            staticnozzleparameters['pcns'], staticnozzleparameters['cstar'],
                            staticnozzleparameters['throatRadiusCurvature'], staticnozzleparameters['at'],
                            gassidearea, Twglist[ind], staticnozzleparameters['tcns'],
                            thrustchamber.machInterpolator(x), params['gamma'])
        hglist[ind] = 1/(1/hglist[ind]+1/staticnozzleparameters['hcoating'])
        qdotlist[ind] = hglist[ind] * (Trlist[ind] - Twglist[ind])
        # This works since its the same method used in qdotdiff function, its just a one iteration approx
        Twclist[ind] = Twglist[ind] - qdotlist[ind] / (matprops['kw']((Twglist[ind])) / twlist[ind])
        Twclist[ind] = Twglist[ind] - qdotlist[ind] / (matprops['kw']((Twglist[ind] + Twclist[ind]) / 2) / twlist[ind])
        hclist[ind] = hc(Twclist[ind], Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                         staticnozzleparameters['pObj'], staticnozzleparameters['pObjWater'],
                         staticnozzleparameters['ethpercent'], coolingfactorlist[ind])
        Tc = Tc + qdotlist[ind] * salist[ind] / (params['mdot_fuel'] * heatCapacity(Tc, staticnozzleparameters))
        rholist[ind] = rho(Tc, coolantpressure, staticnozzleparameters)
        viscositylist[ind] = viscosity(Tc, coolantpressure, staticnozzleparameters['pObj'],
                                       staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])
        Relist[ind] = reynolds(Tc,coolantpressure,vlist[ind],hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                                       staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])

    return Twglist, hglist, qdotlist, Twclist, hclist, Tclist, \
           coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist

