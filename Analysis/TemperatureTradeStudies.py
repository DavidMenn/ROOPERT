import sys
sys.path.insert(1,"./")
import numpy as np
import math
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
from matplotlib import cm
import Analysis.FirstOrderCalcs as FAC
import Components.ThrustChamber as ThrustChamber
import Components.CoolingSystem as CS
import Components.StructuralApproximation as SA
import Toolbox.CADAssistant as CAD
from scipy.optimize import minimize_scalar
import scipy.optimize
import Main
import os

"make program that generates valid corss sections"


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
    heatingfactorlist = np.ones(xlist.size)*.3 # .6 is from cfd last year, i think its bs but whatever
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




#find rli
#find rlist
# preset geometry
# run temps
# run structures
# return FS's



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
        'lstar': 1.24,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa,
        'impulseguess' :  495555.24828424345,
        'rc' : .11} #623919}
configtitle = "NitroEthanol75 at 300 psi Optimize_Wall_Thickness"
params = Main.main(args, configtitle, output=False)

params['thetac'] = (40*math.pi/180)
#conevol = math.pi*params['rc']**3*math.tan(params['thetac'])/3 - math.pi*params['rt']**3*math.tan(params['thetac'])/3
volfunc = lambda lc : math.pi*params['rc']**2*lc  +\
     math.pi*params['rc']**3/math.tan(params['thetac'])/3 -\
         math.pi*params['rt']**3/math.tan(params['thetac'])/3
lstarminimizer = lambda lc : volfunc(lc)/(params['rt']**2*math.pi) - params['lstar']
result = scipy.optimize.root(lstarminimizer, .05, args=(), method='hybr', jac=None, tol=None, callback=None, options=None)
params['lc']=result['x'][0]
xlist = np.linspace(0, params['lc'] + (params['rc'] - params['rt']) / math.tan(params['thetac']) + params['ln_conical'], 25)
       
rlist,xlist = RListGenerator.paraRlist(xlist, params['lc'], params['rc'],
                                  params['lc'] + (params['rc'] - params['rt'])/(math.tan(params['thetac'])),
                                  params['rt'],
                                  params['lc'] + (params['rc'] - params['rt'])/(math.tan(params['thetac'])) + params['ln_conical'],
                                  params['re'], params['lc']*1.8, 2*.0254, 2*.0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
# xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar

TC = ThrustChamber.ThrustChamber(rlist, xlist)
print(TC.rt, TC.xt, TC.xns)

machlist, preslist, templist = TC.flowSimple(params)




# in general, we will do the ith - row is the ith configuration, and the jth column is the jth axial step
numconfigs = 5
xlistflippedmat = np.zeros((numconfigs,len(xlist)))
rlistflippedmat = np.zeros((numconfigs,len(xlist)))
chlistmat  = np.zeros((numconfigs,len(xlist)))
twlistmat  = np.zeros((numconfigs,len(xlist)))
nlistmat  = np.zeros((numconfigs,len(xlist)))
ewlistmat  = np.zeros((numconfigs,len(xlist)))
helicitylistmat  = np.zeros((numconfigs,len(xlist)))
chanelToLandRatiomat = np.zeros(numconfigs)

# YOUR CODE HERE
configsnames = ["trc=.03","trc=.035","trc=.04","trc=.045","trc=.05" ]
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.003
    twlistmat[i,:]  = (rlistflippedmat[0,:]/TC.rt)*.001 
    nlistmat[i,:]  = np.ones(len(xlist))*80
    ewlistmat[i,:]  = chlistmat[0,:]*2
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = (rlistflippedmat[0,:]**math.log(2,params['rc']/params['rt'])/TC.rt**math.log(2,params['rc']/params['rt']))*45*math.pi/180
    chanelToLandRatiomat[i] = 2
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
# this is just to test ideas till I get something I like
"""
configsnames=['helicity**math.log(2,params[rc]/params[rt])',"hel **1.5"]
xlistflippedmat[0,:] = np.flip(xlist)
rlistflippedmat[0,:]  = np.flip(rlist)
chlistmat[0,:]  = (TC.rt/rlistflippedmat[0,:])**.5*.003 
twlistmat[0,:]  = (rlistflippedmat[0,:]/TC.rt)*.001 
nlistmat[0,:]  = np.ones(len(xlist))*80
ewlistmat[0,:]  = np.ones(len(xlist))*.008
#HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
helicitylistmat[0,:]  = (rlistflippedmat[0,:]**math.log(2,params['rc']/params['rt'])/TC.rt**math.log(2,params['rc']/params['rt']))*45*math.pi/180
chanelToLandRatiomat[0] = 2
for index in range(0,np.size(helicitylistmat[0,:] )):
    if helicitylistmat[0,:] [index]>math.pi/2:
        helicitylistmat[0,:] [index] = math.pi/2

xlistflippedmat[1,:] = np.flip(xlist)
rlistflippedmat[1,:]  = np.flip(rlist)
chlistmat[1,:] = (TC.rt/rlistflippedmat[0,:])**.5*.003 
twlistmat[1,:] = (rlistflippedmat[0,:]/TC.rt)*.001 
nlistmat[1,:] = np.ones(len(xlist))*80
ewlistmat[1,:] = np.ones(len(xlist))*.005
#HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
helicitylistmat[1,:]  = (rlistflippedmat[0,:]**1.5/TC.rt**1.5)*45*math.pi/180
chanelToLandRatiomat[1] = 2
for index in range(0,np.size(helicitylistmat[1,:])):
    if helicitylistmat[1,:][index]>math.pi/2:
        helicitylistmat[1,:][index] = math.pi/2
"""
"""
configsnames = ['helicity = 30 deg','helicity = 45 deg','helicity = 60 deg','helicity = 75 deg','helicity = 90 deg']
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.009
    twlistmat[i,:]  = np.ones(len(xlist))*.001
    nlistmat[i,:]  = np.ones(len(xlist))*40
    ewlistmat[i,:]  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = np.ones(len(xlist))*(30+15*i)*math.pi/180
    chanelToLandRatiomat[i] = 2
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
"""
"""
configsnames = ['helicity = 30 deg','helicity = 45 deg','helicity = 60 deg','helicity = 75 deg','helicity = 90 deg']
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.009
    twlistmat[i,:]  = np.ones(len(xlist))*.001
    nlistmat[i,:]  = np.ones(len(xlist))*40
    ewlistmat[i,:]  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = np.ones(len(xlist))*(30+15*i)*math.pi/180
    chanelToLandRatiomat[i] = 2
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
"""
"""
 # this is to compare n
configsnames = ['n = 40','n = 80','n = 160','n=320','n=640']
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.009
    twlistmat[i,:]  = np.ones(len(xlist))*.001
    nlistmat[i,:]  = np.ones(len(xlist))*10*2**(i+2)
    ewlistmat[i,:]  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = np.ones(len(xlist))*math.pi/2
    chanelToLandRatiomat[i] = 2
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
"""
"""
 # this is to compare chaneltoland
configsnames = ['C2L = 1','C2L = 2','C2L = 3','C2L = 4','C2L = 5']
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.009
    twlistmat[i,:]  = np.ones(len(xlist))*.001
    nlistmat[i,:]  = np.ones(len(xlist))*40
    ewlistmat[i,:]  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = np.ones(len(xlist))*math.pi/2
    chanelToLandRatiomat[i] = 1+i
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
"""
"""
 # this is to compare ch
configsnames = ['ch = 4 mm','ch = 6 mm','ch = 8 mm','ch = 10 mm','ch = 12 mm']
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.004+i*.002 
    twlistmat[i,:]  = np.ones(len(xlist))*.001
    nlistmat[i,:]  = np.ones(len(xlist))*40
    ewlistmat[i,:]  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = np.ones(len(xlist))*math.pi/2
    chanelToLandRatiomat[i] = 2
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
"""
"""
 # this is to compare wall thicknesses
configsnames = ['tw = .5 mm','tw = .75 mm','tw = 1 mm','tw = 1.25 mm','tw = 1.6 mm']
for i in range(0,numconfigs):
    xlistflippedmat[i,:] = np.flip(xlist)
    rlistflippedmat[i,:]  = np.flip(rlist)
    chlistmat[i,:]  = np.ones(len(xlist))*.003 
    twlistmat[i,:]  = np.ones(len(xlist))*(.0005+i*.00025)
    nlistmat[i,:]  = np.ones(len(xlist))*40
    ewlistmat[i,:]  = np.ones(len(xlist))*.005
    #HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
    helicitylistmat[i,:]  = np.ones(len(xlist))*math.pi/2
    chanelToLandRatiomat[i] = 2
    for index in range(0,np.size(helicitylistmat[0,:] )):
        if helicitylistmat[i,:] [index]>math.pi/2:
            helicitylistmat[i,:] [index] = math.pi/2
"""
"""
 # this is to compare helicities
xlistflippedmat[0,:] = np.flip(xlist)
rlistflippedmat[0,:]  = np.flip(rlist)
chlistmat[0,:]  = (TC.rt/rlistflippedmat[0,:])**.5*.003 
twlistmat[0,:]  = (rlistflippedmat[0,:]/TC.rt)*.001 
nlistmat[0,:]  = np.ones(len(xlist))*40
ewlistmat[0,:]  = np.ones(len(xlist))*.005
#HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
helicitylistmat[0,:]  = (rlistflippedmat[0,:]**1.5/TC.rt**1.5)*30*math.pi/180
chanelToLandRatiomat[0] = 2
for index in range(0,np.size(helicitylistmat[0,:] )):
    if helicitylistmat[0,:] [index]>math.pi/2:
        helicitylistmat[0,:] [index] = math.pi/2

xlistflippedmat[1,:] = np.flip(xlist)
rlistflippedmat[1,:]  = np.flip(rlist)
chlistmat[1,:] = (TC.rt/rlistflippedmat[1,:])**.5*.003 
twlistmat[1,:] = (rlistflippedmat[0,:]/TC.rt)*.001 
nlistmat[1,:] = np.ones(len(xlist))*40
ewlistmat[1,:] = np.ones(len(xlist))*.005
#HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
helicitylistmat[1,:] = np.ones(len(xlist))*math.pi/2
chanelToLandRatiomat[1] = 2
for index in range(0,np.size(helicitylistmat[1,:])):
    if helicitylistmat[1,:][index]>math.pi/2:
        helicitylistmat[1,:][index] = math.pi/2
"""
vlistflippedmat = np.zeros((numconfigs,len(xlist)))
hydraulicdiamlistmat =  np.zeros((numconfigs,len(xlist)))
cwlistmat = np.zeros((numconfigs,len(xlist)))
Twglistmat = np.zeros((numconfigs,len(xlist)))
hglistmat = np.zeros((numconfigs,len(xlist)))
qdotlistmat = np.zeros((numconfigs,len(xlist)))
Twclistmat = np.zeros((numconfigs,len(xlist)))
hclistmat = np.zeros((numconfigs,len(xlist)))
Tclistmat = np.zeros((numconfigs,len(xlist)))
coolantpressurelistmat = np.zeros((numconfigs,len(xlist)))
fincoolingfactorlistmat = np.zeros((numconfigs,len(xlist)))
rholistmat = np.zeros((numconfigs,len(xlist)))
viscositylistmat = np.zeros((numconfigs,len(xlist)))
Relistmat = np.zeros((numconfigs,len(xlist)))
FOSlistmat = np.zeros((numconfigs,len(xlist)))
for i in range(0,numconfigs):
    params['throat_radius_curvature'] = .03+i*.005
    print(f"trc = {params['throat_radius_curvature']}")
    alistflipped, xlist, vlistflippedmat[i,:],\
        hydraulicdiamlistmat[i,:], salistflipped, dxlist, fincoolingfactorfunc, cwlistmat[i,:],\
            Twglistmat[i,:], hglistmat[i,:], qdotlistmat[i,:], Twclistmat[i,:], hclistmat[i,:], Tclistmat[i,:], coolantpressurelistmat[i,:], qdotlist, fincoolingfactorlistmat[i,:], rholistmat[i,:], viscositylistmat[i,:], Relistmat[i,:], FOSlistmat[i,:] = RunCoolingSystem(chlistmat[i,:],
    twlistmat[i,:],
    nlistmat[i,:],
    helicitylistmat[i,:],
    params,
    xlist,
    rlist,
    chanelToLandRatiomat[i], 
    TC )

plt.figure()
for i in range(0,numconfigs):
    plt.plot(xlistflippedmat[i,:],FOSlistmat[i,:],label = f"{configsnames[i]}, dT = {int(np.max(Twglistmat[i,:]-Twclistmat[i,:]))} K")
plt.xlabel("TC position, meters")
plt.ylabel("Factor of Safety")
plt.legend()
plt.title("Factor of Safeties")

plt.figure()
for i in range(0,numconfigs):
    plt.plot(xlistflippedmat[i,:],Twglistmat[i,:],label = configsnames[i])
plt.xlabel("TC position, meters")
plt.ylabel("Wall Gas Temp [K]")
plt.legend()
dxtrunc=(xlistflippedmat[0,0]-xlistflippedmat[0,1])
#plt.title(f"Gas Side Wall Temps, steps ={len(xlistflippedmat[0,:])}, dx = {'%.5f'%dxtrunc}")
plt.title(f"Gas Side Wall Temps")

plt.figure()
for i in range(0,numconfigs):
    plt.plot(xlistflippedmat[i,:],Twclistmat[i,:],label = configsnames[i])
plt.xlabel("TC position, meters")
plt.ylabel("Coolant SideWall Gas Temp [K]")
plt.legend()
plt.title("Coolant Side Wall Temps")

plt.figure()
for i in range(0,numconfigs):
    plt.plot(xlistflippedmat[i,:],fincoolingfactorlistmat[i,:],label = configsnames[i])
plt.xlabel("TC position, meters")
plt.ylabel('Fin Factor (Corrected coolant side area/gas side area)')
plt.legend()
plt.title("Fin Factors")

plt.figure()
for i in range(0,numconfigs):
    plt.plot(xlistflippedmat[i,:],Tclistmat[i,:],label = configsnames[i])
plt.xlabel("TC position, meters")
plt.ylabel("Coolant Temp [K]")
plt.legend()
plt.title("Coolant Temps")

plt.figure()
for i in range(0,numconfigs):
    plt.plot(xlistflippedmat[i,:],coolantpressurelistmat[i,:]/const.psiToPa,label = configsnames[i])
plt.xlabel("TC position, meters")
plt.ylabel("Coolant pressure [psi]")
plt.legend()
plt.title("Coolant pressures")

plt.figure()
for i in range(0,numconfigs):
    plt.plot(np.hstack((xlistflippedmat[i,:],xlist)),np.hstack((rlistflippedmat[i,:],-rlist)),'k',label = configsnames[i])
plt.xlabel("TC position [m]")
plt.ylabel("Radius [m]")
plt.title("Thrust Chamber Shape")


for i in range(0,numconfigs):
    plt.figure()
    plt.plot(xlistflippedmat[i,:],chlistmat[i,:]*1000,'r',label="Chanel Height [mm]")
    plt.plot(xlistflippedmat[i,:],cwlistmat[i,:]*1000,'b',label="Chanel Width [mm]")
    plt.plot(xlistflippedmat[i,:],twlistmat[i,:]*1000,'k',label="Wall Thickness [mm]")
    plt.plot(xlistflippedmat[i,:],hydraulicdiamlistmat[i,:]*1000,'g',label="Hydraulic Diam [mm]")
    plt.plot(xlistflippedmat[i,:],helicitylistmat[i,:]*180/math.pi/10,'m',label="helicity [10's of degrees]")
    plt.plot(xlistflippedmat[i,:],vlistflippedmat[i,:],'c',label="Coolant Velocity [m/s]")
    plt.legend()
    plt.xlabel("TC position, meters")
    plt.ylabel("Thicknesses, [mm]")
    plt.title(f"Geometry for {configsnames[i]}, n = {int(nlistmat[i,1])}, Chanel to Land = {chanelToLandRatiomat[i]}")

plt.show()
"""
title="ChamberTemps"
fig, axs = plt.subplots(3,3)
fig.suptitle(title)

axs[0,1].plot(xlistflipped,hglist , 'g')  # row=0, column=0
axs[1,1].plot(xlistflipped,hclist , 'r')# row=1, column=0
axs[2,1].plot(np.hstack((xlistflipped,xlist)),np.hstack((rlistflipped,-rlist)) , 'k') # row=0, column=0


axs[0,1].set_title('hglist')
axs[1,1].set_title('hclist')
axs[2,1].set_title('Thrust Chamber Shape')

axs[0,0].plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp")
axs[0,0].plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp") # row=0, column=0
axs[0,0].plot(xlistflipped,Tclist , 'b', label="Coolant Temp") #
axs[1,0].plot(xlistflipped,Tclist , 'r')# row=1, column=0
axs[2,0].plot(xlistflipped,fincoolingfactorlist/(2*math.pi*rlist/nlist), 'r')# row=1, column=0

axs[0,0].set_title('Twg')
axs[1,0].set_title('Tc')
axs[2,0].set_title('Fin Factor (Corrected coolant side area/gas side area)')
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
print(f"max twg = {np.max(Twglist)} in kelvin, {np.max(Twglist)*const.degKtoR} in Rankine (freedom)\n max Twc ="
        f" {np.max(Twclist)} in kelvin, {np.max(Twclist)*const.degKtoR} in Rankine (freedom)")

plt.figure()
plt.plot(xlistflipped,chlist,'r',label="Chanel Height [mm]")
plt.plot(xlistflipped,cwlist,'b',label="Chanel Widht [mm]")
plt.plot(xlistflipped,twlist,'k',label="Wall Thickness [mm]")
plt.plot(xlistflipped,hydraulicdiamlist,'g',label="Hydraulic Diam [mm]")
plt.legend()
plt.xlabel("TC position, meters")
plt.ylabel("thicknesses, meters")
plt.title("Cooling Passage Geometry")


plt.figure()
plt.plot(xlistflipped,FOSlist,'k')
plt.xlabel("TC position, meters")
plt.ylabel("Factor of Safety")
plt.title("Factor of Safeties")

plt.figure()
plt.plot(xlistflipped,viscositylist/rholist,'k')
plt.xlabel("TC position, meters")
plt.ylabel("Kinematic viscosity (Pa*s*m^3/kg)")
plt.title("Kinematic Viscosity")

plt.figure()
plt.plot(Tclist,viscositylist/rholist,'k')
plt.xlabel("Coolant Temp, Kelvin")
plt.ylabel("Kinematic viscosity (Pa*s*m^3/kg)")
plt.title(f"Kinematic Viscosity vs Temperature @P_coolant={math.floor(coolantpressurelist[0])}")

plt.show()
"""