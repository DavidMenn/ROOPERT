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

class StructuralAnalysis(object): #This object should be used to calculate stersses and factors of safety given geometry and temperatures

    def __init__(self, rlist, xlist, nlist, chlist, cwlist, twlist, material): #Initializes geometry and useful values for rectangular changels
        self.rlist = rlist
        self.xlist = xlist
        self.nlist = nlist
        self.chlist = chlist
        self.cwlist=cwlist    
        self.cwlist=cwlist    
        self.twlist=twlist

        self.rclist = np.zeros(rlist.size) #Radius of curvature is defined such that a negative readius of curvature is convex, positive is concave (like the throat)
        ind = 1
        dr2dx2 = (rlist[ind+1]-2*rlist[ind]+rlist[ind-1])/(((xlist[ind-1]-xlist[ind+1])/2)**2)
        if dr2dx2 == 0:
            self.rclist[ind]=0
        else:
            drdx = ((rlist[ind]-rlist[ind-1])/(xlist[ind]-xlist[ind-1])+(rlist[ind]-rlist[ind-+1])/(xlist[ind]-xlist[ind+1]))/2
            self.rclist[ind]=((1+drdx**2)**(3/2))/(dr2dx2)
        self.rclist[0]=self.rclist[1]
        while ind<xlist.size-1:
            dr2dx2 = (rlist[ind+1]-2*rlist[ind]+rlist[ind-1])/(((xlist[ind-1]-xlist[ind+1])/2)**2)
            if dr2dx2 == 0:
                self.rclist[ind]=0
            else:
                drdx = ((rlist[ind]-rlist[ind-1])/(xlist[ind]-xlist[ind-1])+(rlist[ind]-rlist[ind-+1])/(xlist[ind]-xlist[ind+1]))/2
                self.rclist[ind]=((1+drdx**2)**(3/2))/(dr2dx2)
            ind=ind+1
        self.rclist[-1]=self.rclist[-2]

        #CURRENTLY HAVE THIS AS FIXED, NEED TO BE ABLE TO CHANGE THIS WITH PRELOADED MATPROPS
        self.EInterpolator = interpolate.interp1d(((np.hstack((70,np.arange(100, 2001, 100), np.array([10000])))-32)*5/9)+273.15,
                                              np.array([29,28.8,28.4,28,27.6,27.1,26.7,26.3,25.8,25.3,24.8,24.2,23.7,23.0,22.3,21.3,20.2,18.8,17.4,15.9,14.3,0])*(10**6)*const.psiToPa, kind='linear') #https://www.engineersedge.com/materials/inconel_718_modulus_of_elasticity_table_13562.htm
        self.PoissonInterpolator = interpolate.interp1d(((np.hstack((70,np.arange(100, 2001, 100), np.array([10000])))-32)*5/9)+273.15,
                                              np.array([.294,.291,.288,.280,.275,.272,.273,.271,.272,.271,.276,.283,.292,.306,.321,.331,.334,.341,.366,.402,.402,.402]), kind='linear') #https://www.engineersedge.com/materials/inconel_718_modulus_of_elasticity_table_13562.htm
        #self.SyInterpolator = interpolate.interp1d(np.array([25,425,650,870,871,950,10000]),
        #                                      np.array([1200,1050,1000,300,290,100,1])*10**6, kind='linear') # pg 734 metal additive manufacturing, inconel 718
        self.SyInterpolator = interpolate.interp1d(np.array([25,425,650,870,871,10000])+273,
                                              np.array([1200,1050,1000,300,10,1])*10**6, kind='linear') # pg 734 metal additive manufacturing, inconel 718
        
        self.AlphaInterpolator = interpolate.interp1d(np.hstack((np.arange(0, 1001, 100), np.array([10000])))+273,
                                              np.array([12.6,12.6, 13.9, 14.2, 14.5, 14.8, 15.1, 15.6, 16.4, 17.5, 17.8,17.8])*10**(-6), kind='linear') # pg 781 metal additive manufacturing, inconel 718

    def FOS(self, Twglist, Twclist, coolantpressurelist, preslist): #Calculates stresses using just sigma_tangential from heister page 207
        dplist = coolantpressurelist-preslist
        avgtemplist = (Twglist+Twclist)/2
        deltaTlist = Twglist-Twclist
        Elist = self.EInterpolator(avgtemplist)
        poissonlist = self.PoissonInterpolator(avgtemplist)
        Sylist = self.SyInterpolator(avgtemplist)
        Alphalist = self.AlphaInterpolator(avgtemplist)
        FUDGEFACTORFORTHERMALSTRESS = 1#.02 #pretty much disregarding it rn lol
        sigmatangentiallist = dplist/2*((self.cwlist/self.twlist)**2)   +  FUDGEFACTORFORTHERMALSTRESS*Elist*Alphalist*deltaTlist/(2*(1-poissonlist))
        return Sylist/sigmatangentiallist


    


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



def steadyStateTemperatures(wallmaterial, thrustchamber, params, salist, nlist, coolingfactorlist,
                            heatingfactorlist, xlist, vlist, initialcoolanttemp, initialcoolantpressure, twlist,
                            hydraulicdiamlist, rgaslist = None, fincoolingfactorfunc = None, dxlist = None, helicitylist = None):
    if fincoolingfactorfunc is None: # This is so old code with coaxial shell still runs, ideally you use Qdot so you can use thermal resistances!
        useqdot=True
        if type(nlist) is int:
            nlist = nlist*np.ones(xlist.size)
    else:
        useqdot=False #this is new!

    if helicitylist is None:
        helicitylist = np.ones(np.size(xlist))*math.pi/2

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
        kwInterpolator = interpolate.interp1d(np.hstack((-1000,25,np.arange(100, 1001, 100) + 273.5, np.array([10000]))),
                                              np.array([0,12.9, 13.9, 16.1, 18.2, 20, 22.1, 24.7, 31.4, 30.1, 31.5, 36.7,
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
        'roughness' : 320*10**-6,
        'nlist' : nlist,
        'helicitylist' : helicitylist,
        'coatingthickness' : .0005,
        'kcoating' : 10.1,
        'hcoating': 14470
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
    if not useqdot:
        Qdotlist = np.zeros(xlist.size)
        fincoolingfactorlist = np.zeros(xlist.size)
    Twclist = np.zeros(xlist.size)
    hclist = np.zeros(xlist.size)
    Tclist = np.zeros(xlist.size)
    rholist = np.zeros(xlist.size)
    viscositylist = np.zeros(xlist.size)
    Relist = np.zeros(xlist.size)
    
    # THIS IS JUST TO output CF plot
    """
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
    if useqdot: # THIS SHOULD NEVER BE USED!! OLD VERSION
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
            if ind%50 == 0:
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
            Tc = Tc + qdotlist[ind] * salist[ind] / (params['mdot_fuel']/nlist[ind] * heatCapacity(Tc, staticnozzleparameters))
            rholist[ind] = rho(Tc, coolantpressure, staticnozzleparameters)
            viscositylist[ind]= viscosity(Tc,coolantpressure,  staticnozzleparameters['pObj'],staticnozzleparameters['pObjWater'],staticnozzleparameters['ethpercent'])
            Relist[ind] = reynolds(Tc,coolantpressure,vlist[ind],hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                                        staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])
        return Twglist, hglist, qdotlist, Twclist, hclist, Tclist, \
           coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist
    else:
        for ind in np.arange(0, xlist.size):
            Tclist[ind] = Tc  # set cooland at current station to coolant temp
            dx = 0
            if ind == 0:
                dx = xlist[ind] - xlist[ind + 1]
                coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                        staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                        dx) / \
                                hydraulicdiamlist[ind] * (
                                            .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)
            else:
                dx = xlist[ind - 1] - xlist[ind]
                coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                        staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                        dx) / \
                                hydraulicdiamlist[ind] * (
                                        .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)


            coolantpressurelist[ind] = coolantpressure
            x = xlist[ind]
            if ind%10 == 0:
                print(x)
            Tri = Trlist[ind]
            gassidearea = thrustchamber.areaInterpolator(x)
            try:
                if ind == 0:
                    if abs(rgaslist[ind]-params['re'])>abs(rgaslist[-1]-params['re']):
                        rgaslist=np.flip(rgaslist) #you didnt do your list the right way idiot
            except:
                pass #no params['re']
            rgas = rgaslist[ind]
            fincoolingfactorfunc_atstation = lambda hc, kw : fincoolingfactorfunc(hc,kw,ind)

            Twglist[ind] = QdotdiffMinimizer(staticnozzleparameters, gassidearea, thrustchamber.machInterpolator(x),
                                            params['gamma'],
                                            Tri, Tc, twlist[ind], coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                            coolingfactorlist[ind],heatingfactorlist[ind], rgas, nlist[ind], fincoolingfactorfunc_atstation,
                                            helicity = staticnozzleparameters['helicitylist'][ind])
            hglist[ind] = heatingfactorlist[ind]*Bartz(staticnozzleparameters['throatdiameter'], staticnozzleparameters['viscosityns'],
                                staticnozzleparameters['prandtlns'], staticnozzleparameters['cpns'],
                                staticnozzleparameters['pcns'], staticnozzleparameters['cstar'],
                                staticnozzleparameters['throatRadiusCurvature'], staticnozzleparameters['at'],
                                gassidearea, Twglist[ind], staticnozzleparameters['tcns'],
                                thrustchamber.machInterpolator(x), params['gamma'])
            Qdotlist[ind] = hglist[ind] * (Trlist[ind] - Twglist[ind])*2*math.pi*rgas/nlist[ind]
            qdotlist[ind] = hglist[ind] * (Trlist[ind] - Twglist[ind])
            # This works since its the same method used in qdotdiff function, its just a one iteration approx
            Twclist[ind] = Twglist[ind] - Qdotlist[ind]/(2*math.pi*rgas/nlist[ind]) / (matprops['kw']((Twglist[ind])) / twlist[ind])
            Twclist[ind] = Twglist[ind] - Qdotlist[ind]/(2*math.pi*rgas/nlist[ind]) / (matprops['kw']((Twglist[ind]+Twclist[ind])/2) / twlist[ind])
            hclist[ind] = hc(Twclist[ind], Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                            staticnozzleparameters['pObj'], staticnozzleparameters['pObjWater'],
                            staticnozzleparameters['ethpercent'], coolingfactorlist[ind])
            fincoolingfactorlist[ind] = dx*fincoolingfactorfunc_atstation(hclist[ind],(matprops['kw']((Twglist[ind]+Twclist[ind])/2)))\
                                        /math.sin(staticnozzleparameters['helicitylist'][ind])/(gassidearea/nlist[ind])
            #This is the ratio of corrected coolant side area over gas side area servicied by coolant passage
            #This would be (2*ch*finefficiency + cw)*(dx/sin(helicity)) / (gassidearea/numchanels)
            Tc = Tc + Qdotlist[ind] * dxlist[ind] / (params['mdot_fuel']/nlist[ind] * heatCapacity(Tc, staticnozzleparameters))
            rholist[ind] = rho(Tc, coolantpressure, staticnozzleparameters)
            viscositylist[ind]= viscosity(Tc,coolantpressure,  staticnozzleparameters['pObj'],staticnozzleparameters['pObjWater'],staticnozzleparameters['ethpercent'])
            Relist[ind] = reynolds(Tc,coolantpressure,vlist[ind],hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                                        staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])
            
        return Twglist, hglist, Qdotlist, Twclist, hclist, Tclist, \
           coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Trlist # lol this is bandaid should be Relist

def ablative():
    return "set this up"


def transientTemperature(wallmaterial, thrustchamber, params, salist, nlist, coolingfactorlist,
                            heatingfactorlist, xlist, vlist, initialcoolanttemp, initialcoolantpressure, twlist,
                            hydraulicdiamlist, rgaslist = None, fincoolingfactorfunc = None, dxlist = None, helicitylist = None):
    if fincoolingfactorfunc is None: # This is so old code with coaxial shell still runs, ideally you use Qdot so you can use thermal resistances!
        useqdot=True
        if type(nlist) is int:
            nlist = nlist*np.ones(xlist.size)
    else:
        useqdot=False #this is new!

    if helicitylist is None:
        helicitylist = np.ones(np.size(xlist))*math.pi/2

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
        kwInterpolator = interpolate.interp1d(np.hstack((-1000,25,np.arange(100, 1001, 100) + 273.5, np.array([10000]))),
                                              np.array([0,12.9, 13.9, 16.1, 18.2, 20, 22.1, 24.7, 31.4, 30.1, 31.5, 36.7,
                                                        36.7]), kind='linear')
        matprops = {'kw': kwInterpolator}  # copper = 398, mild steel = 51??
        matprops['rho'] = 8220 #kg/m3 for inconcel 718

        cpInterpolator = interpolate.interp1d(np.hstack((25,np.arange(100, 1001, 100) + 273.5, np.array([10000]))),
                                             1000* np.array([.537,.543,.577,.596,.608,.628,.665,.722,.756,.785,.874,.874]), kind='linear')
        matprops['cp'] = cpInterpolator

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
        'roughness' : 320*10**-6,
        'nlist' : nlist,
        'helicitylist' : helicitylist,
        'coatingthickness' : .0005,
        'kcoating' : 10.1,
        'hcoating': 144700,
        'rho': matprops['rho'],
        'cp' : cpInterpolator
    }

    Vollist = 2*math.pi*rgaslist/nlist*dxlist*twlist #volume of the metal that is heating up per chanel


    Trlist = np.zeros(xlist.size)
    for ind in np.arange(0, xlist.size):
        x = xlist[ind]
        Trlist[ind] = recoveryTemp(thrustchamber.tempInterpolator(xlist[ind]), params['gamma'],
                                   thrustchamber.machInterpolator(xlist[ind]), Pr=params['pr_throat'])

    Tc = initialcoolanttemp
    coolantpressure = initialcoolantpressure
    coolantpressurelist = np.zeros(xlist.size)
    Twglist = np.ones(xlist.size)*initialcoolanttemp #initializing the gas side wall temp at room temp
    hglist = np.zeros(xlist.size)
    qdotlist = np.zeros(xlist.size)
    if not useqdot:
        Qdotlist = np.zeros(xlist.size)
        fincoolingfactorlist = np.zeros(xlist.size)
    Twclist = np.ones(xlist.size)*initialcoolanttemp #initializing the coolant side wall temp at room temp
    hclist = np.zeros(xlist.size)
    Tclist = np.zeros(xlist.size)
    rholist = np.zeros(xlist.size)
    viscositylist = np.zeros(xlist.size)
    Relist = np.zeros(xlist.size)
    
    maxtime = 100# seconds Max time
    time = 0
    dt=.1

    twgprevious = -100 
    error_in_steadystate = .1 #measured in change in degrees/second

    while abs(Twglist[round(xlist.size/2)]-twgprevious)/dt > error_in_steadystate and time < maxtime:
        twgprevious = Twglist[round(xlist.size/2)]#pick a random index to check if it gets to steady state

    

    

        if useqdot: # THIS SHOULD NEVER BE USED!! OLD VERSION
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
                if ind%50 == 0:
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
                Tc = Tc + qdotlist[ind] * salist[ind] / (params['mdot_fuel']/nlist[ind] * heatCapacity(Tc, staticnozzleparameters))
                rholist[ind] = rho(Tc, coolantpressure, staticnozzleparameters)
                viscositylist[ind]= viscosity(Tc,coolantpressure,  staticnozzleparameters['pObj'],staticnozzleparameters['pObjWater'],staticnozzleparameters['ethpercent'])
                Relist[ind] = reynolds(Tc,coolantpressure,vlist[ind],hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                                            staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])
            return Twglist, hglist, qdotlist, Twclist, hclist, Tclist, \
            coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist
        else:

            for ind in np.arange(0, xlist.size):
                Tclist[ind] = Tc  # set cooland at current station to coolant temp
                dx = 0
                if ind == 0:
                    dx = xlist[ind] - xlist[ind + 1]
                    coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                            staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                            dx) / \
                                    hydraulicdiamlist[ind] * (
                                                .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)
                else:
                    dx = xlist[ind - 1] - xlist[ind]
                    coolantpressure = coolantpressure - 4 * cf(Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind],
                                                            staticnozzleparameters, roughness=staticnozzleparameters['roughness']) * (
                                            dx) / \
                                    hydraulicdiamlist[ind] * (
                                            .5 * rho(Tc, coolantpressure, staticnozzleparameters) * vlist[ind] ** 2)


                coolantpressurelist[ind] = coolantpressure
                x = xlist[ind]
                
                Tri = Trlist[ind]
                gassidearea = thrustchamber.areaInterpolator(x)
                try:
                    if ind == 0:
                        if abs(rgaslist[ind]-params['re'])>abs(rgaslist[-1]-params['re']):
                            rgaslist=np.flip(rgaslist) #you didnt do your list the right way idiot
                except:
                    pass #no params['re']
                rgas = rgaslist[ind]
                fincoolingfactorfunc_atstation = lambda hc, kw : fincoolingfactorfunc(hc,kw,ind)

                hglist[ind] = heatingfactorlist[ind]*Bartz(staticnozzleparameters['throatdiameter'], staticnozzleparameters['viscosityns'],
                                    staticnozzleparameters['prandtlns'], staticnozzleparameters['cpns'],
                                    staticnozzleparameters['pcns'], staticnozzleparameters['cstar'],
                                    staticnozzleparameters['throatRadiusCurvature'], staticnozzleparameters['at'],
                                    gassidearea, Twglist[ind], staticnozzleparameters['tcns'],
                                    thrustchamber.machInterpolator(x), params['gamma'])

                Qdotlist[ind] = hglist[ind] * (Tri - Twglist[ind])*2*math.pi*rgas/nlist[ind] #This is Qdot per unit length, so it should be multiplied by dx to get actual total heat flux. This overcomplicates solver
                Qdot_throughwall = abs(staticnozzleparameters['kw']((Twglist[ind]+Twclist[ind])/2)/twlist[ind]*(Twglist[ind]-Twclist[ind]))*2*math.pi*rgas/nlist[ind] 
               
                #update twglist by mupltipliyng Qdot by dt, then divide by heat capacity and mass CURRENTLY ASSUMING MASS IS HALF OF TOTAL FOR CHANEL
                Twglist[ind] =Twglist[ind] + dt*(Qdotlist[ind] - Qdot_throughwall)*dxlist[ind] /(Vollist[ind]/2*matprops['rho']* matprops['cp']((Twglist[ind]+Twclist[ind])/2))

                #Check huzel and huang page 98 for graph of effect of helicity on hc
                turningangle = math.pi/2-helicitylist[ind]
                curvature_enhancement_factor =  np.max([1,np.min([1.4/18*turningangle*180/math.pi-1.4/18*20+1, 1.4,
                -1.4/20*turningangle*180/math.pi+1.4/20*80])])#cheesy linear interp
                hci = curvature_enhancement_factor*hc(Twclist[ind], Tc, coolantpressure, vlist[ind], hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                        staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'], coolingfactorlist[ind]) #only reason we guess twci now is to get hci, otherwise we are calculating the thermal resistance directly from twg to tc
                nf = fincoolingfactorfunc_atstation(hci,staticnozzleparameters['kw']((Twglist[ind]+Twclist[ind])/2)) # THIS IS FIN EFFICIENCY * 2 * CH + CW , if you dont want to factor in fin cooling just set this equal to cw
                Rtotal = (twlist[ind] / staticnozzleparameters['kw']((Twglist[ind]+Twclist[ind])/2)       )*2*math.pi*rgas/nlist[ind]+1/(nf*hci)*math.sin(helicitylist[ind]) #missing a factor of 1/dx, to be consistent with Qdotguess being per unit length
                # dividing by the sin of helicity here implies that the fin is acting over dx/sin(helicity), which it is!
                Qdothc = 1/Rtotal * (Twclist[ind]-Tc) #/ math.sin(helicity) #smaller helicity means longer chanel length per dx, resulting in more heat flux into coolant
                
                Twclist[ind] =Twclist[ind] + dt*(Qdot_throughwall - Qdothc)*dxlist[ind] /(Vollist[ind]/2*matprops['rho']* matprops['cp']((Twglist[ind]+Twclist[ind])/2))

                qdotlist[ind] = hglist[ind] * (Trlist[ind] - Twglist[ind])
                # This works since its the same method used in qdotdiff function, its just a one iteration approx
                
                hclist[ind] = hci
                fincoolingfactorlist[ind] = dx*fincoolingfactorfunc_atstation(hclist[ind],(matprops['kw']((Twglist[ind]+Twclist[ind])/2)))\
                                            /math.sin(staticnozzleparameters['helicitylist'][ind])/(gassidearea/nlist[ind])
                #This is the ratio of corrected coolant side area over gas side area servicied by coolant passage
                #This would be (2*ch*finefficiency + cw)*(dx/sin(helicity)) / (gassidearea/numchanels)
                Tc = Tc +  Qdothc * dxlist[ind] / (params['mdot_fuel']/nlist[ind] * heatCapacity(Tc, staticnozzleparameters))
                rholist[ind] = rho(Tc, coolantpressure, staticnozzleparameters)
                viscositylist[ind]= viscosity(Tc,coolantpressure,  staticnozzleparameters['pObj'],staticnozzleparameters['pObjWater'],staticnozzleparameters['ethpercent'])
                Relist[ind] = reynolds(Tc,coolantpressure,vlist[ind],hydraulicdiamlist[ind], staticnozzleparameters['pObj'],
                                            staticnozzleparameters['pObjWater'], staticnozzleparameters['ethpercent'])
                
        
        
        
        if time==0:
            coolantpressuremat =np.reshape( coolantpressurelist,(1,xlist.size))
            Twgmat =np.reshape(  Twglist,(1,xlist.size))
            hgmat = np.reshape( hglist,(1,xlist.size))
            qdotmat =np.reshape( qdotlist,(1,xlist.size)) 
            if not useqdot:
                Qdotmat =np.reshape(Qdotlist,(1,xlist.size))
                fincoolingfactormat =np.reshape( fincoolingfactorlist,(1,xlist.size))
            Twcmat = np.reshape( Twclist,(1,xlist.size)) 
            hcmat =np.reshape( hclist,(1,xlist.size))
            Tcmat = np.reshape( Tclist,(1,xlist.size))
            rhomat = np.reshape(rholist,(1,xlist.size)) 
            viscositymat =  np.reshape(viscositylist,(1,xlist.size)) 
            Remat =  np.reshape(Relist,(1,xlist.size))
        
        else:
            coolantpressuremat = np.concatenate( (coolantpressuremat,np.reshape( coolantpressurelist,(1,xlist.size))))
            Twgmat =  np.concatenate((Twgmat,np.reshape(  Twglist,(1,xlist.size))))
            hgmat = np.concatenate( (hgmat, np.reshape( hglist,(1,xlist.size))))
            qdotmat =np.concatenate((qdotmat, np.reshape(  qdotlist,(1,xlist.size)) ))
            if not useqdot:
                Qdotmat =np.concatenate( (Qdotmat,np.reshape( Qdotlist,(1,xlist.size))))
                fincoolingfactormat = np.concatenate( (fincoolingfactormat,np.reshape( fincoolingfactorlist,(1,xlist.size))))
            Twcmat = np.concatenate((Twcmat , np.reshape(  Twclist,(1,xlist.size)) ))
            hcmat = np.concatenate((hcmat, np.reshape( hclist,(1,xlist.size))))
            Tcmat = np.concatenate( (Tcmat,np.reshape(  Tclist,(1,xlist.size))))
            rhomat =np.concatenate( ( rhomat, np.reshape( rholist,(1,xlist.size)) ))
            viscositymat =np.concatenate((viscositymat,  np.reshape(  viscositylist,(1,xlist.size)) ))
            Remat = np.concatenate( (Remat, np.reshape( Relist,(1,xlist.size))))
        
        
        Tc = initialcoolanttemp
        coolantpressure = initialcoolantpressure
        coolantpressurelist = np.zeros(xlist.size)
        #Twglist = np.zeros(xlist.size)
        hglist = np.zeros(xlist.size)
        qdotlist = np.zeros(xlist.size)
        if not useqdot:
            Qdotlist = np.zeros(xlist.size)
            fincoolingfactorlist = np.zeros(xlist.size)
        #Twclist = np.zeros(xlist.size)
        hclist = np.zeros(xlist.size)
        Tclist = np.zeros(xlist.size)
        rholist = np.zeros(xlist.size)
        viscositylist = np.zeros(xlist.size)
        Relist = np.zeros(xlist.size)

        print(time)
        time += dt

    return Twgmat, hgmat, Qdotmat, Twcmat, hcmat, Tcmat, \
        coolantpressuremat, qdotmat, fincoolingfactormat, rhomat, viscositymat, Remat, time


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
    return scipy.optimize.minimize_scalar(qdotdiff, bounds=(350, Tri), tol=10, method='bounded',
                                          args=(staticnozzleparams, a, mach, gamma,
                                                Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam,
                                                coolingfactor, heatingfactor))['x']  # going to find root of qdotdiff for twg
def QdotdiffMinimizer(staticnozzleparams, a, mach, gamma,
                      Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam, coolingfactor, heatingfactor, rgas, n, fincoolingfactorfunc,
                      helicity = math.pi/2):
    return scipy.optimize.minimize_scalar(Qdotdiff, bounds=(350, Tri), tol=.01, method='bounded',
                                          args=(staticnozzleparams, a, mach, gamma,
                                                Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam,
                                                coolingfactor, heatingfactor, rgas, n, fincoolingfactorfunc, helicity))['x']  # going to find root of qdotdiff for twg


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


def Qdotdiff(Twgi, staticnozzleparams, a, mach, gamma,
             Tri, Tc, tw, coolantpressure, coolantvelocity, hydraulicdiam, coolingfactor, heatingfactor, rgas, n, fincoolingfactorfunc,
             helicity = math.pi/2):
    hgi = heatingfactor*Bartz(staticnozzleparams['throatdiameter'], staticnozzleparams['viscosityns'],
                staticnozzleparams['prandtlns'], staticnozzleparams['cpns'],
                staticnozzleparams['pcns'], staticnozzleparams['cstar'], staticnozzleparams['throatRadiusCurvature'],
                staticnozzleparams['at'],
                a, Twgi, staticnozzleparams['tcns'], mach, gamma)
    hgi = 1/( 1/hgi) #lol
    Qdotguess = hgi * (Tri - Twgi)*2*math.pi*rgas/n #This is Qdot per unit length, so it should be multiplied by dx to get actual total heat flux. This overcomplicates solver
    Twci = Twgi - ((Qdotguess/(2*math.pi*rgas/n) * tw) / staticnozzleparams['kw'](Twgi))# get the initial guess
    kwi = staticnozzleparams['kw']((Twgi+Twci)/2)
    Twci = Twgi - ((Qdotguess/(2*math.pi*rgas/n) * tw) / staticnozzleparams['kw']((Twgi+Twci)/2)) # now use avg wall temp
    
    #Check huzel and huang page 98 for graph of effect of helicity on hc
    turningangle = math.pi/2-helicity
    curvature_enhancement_factor =  np.max([1,np.min([1.4/18*turningangle*180/math.pi-1.4/18*20+1, 1.4,
     -1.4/20*turningangle*180/math.pi+1.4/20*80])])#cheesy linear interp
    hci = curvature_enhancement_factor*hc(Twci, Tc, coolantpressure, coolantvelocity, hydraulicdiam, staticnozzleparams['pObj'],
             staticnozzleparams['pObjWater'], staticnozzleparams['ethpercent'], coolingfactor) #only reason we guess twci now is to get hci, otherwise we are calculating the thermal resistance directly from twg to tc
    nf = fincoolingfactorfunc(hci,kwi) # THIS IS FIN EFFICIENCY * 2 * CH + CW , if you dont want to factor in fin cooling just set this equal to cw
    Rtotal = (tw / (kwi)       )*2*math.pi*rgas/n+1/(nf*hci)*math.sin(helicity) #missing a factor of 1/dx, to be consistent with Qdotguess being per unit length
    # dividing by the sin of helicity here implies that the fin is acting over dx/sin(helicity), which it is!
    Qdothc = 1/Rtotal * (Twci-Tc) #/ math.sin(helicity) #smaller helicity means longer chanel length per dx, resulting in more heat flux into coolant
    #print(f"{[Twgi, hgi, hci, nf, Qdotguess, Qdothc, abs(Qdothc - Qdotguess)]},")
    return abs(Qdothc - Qdotguess)


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
    # FIGURE OUT THESE COEFFICIENTS PETE, PAGE 197 IN HEISTER eider - tate equation
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
            -2 * math.log(roughness / diameter / 3.7 + 2.51 / Re / math.sqrt(cf * 4 + .00001), 10) - 1 / math.sqrt(
                (cf + .00001) * 4))
def viscosity(    temp, pres,  pObj, pObjWater, ethpercent):
    if ethpercent is None:
        return .1 * pObj.Visc_compressed(temp * const.degKtoR, #Converting from poise to pascal-seconds
                                              pres / const.psiToPa)
    else:

        viscosityeth = .1 * pObj.Visc_compressed(temp * const.degKtoR,
                                                 pres / const.psiToPa)
        viscositywater = .1 * pObjWater.Visc_compressed(temp * const.degKtoR,
                                                        pres / const.psiToPa)

        return  ethpercent * viscosityeth + (1 - ethpercent) * viscositywater

