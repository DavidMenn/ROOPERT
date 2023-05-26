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
    heatingfactorlist = np.ones(xlist.size)*.6 #.3 # .6 is from cfd last year, i think its bs but whatever
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

    Twglist, hglist, Qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist = CS.steadyStateTemperatures(
        None, TC, params, salistflipped, nlist, coolingfactorlist,
        heatingfactorlist, xlist, vlistflipped, 293, params['pc'] + params['pc']*.2 + 50*const.psiToPa, twlistflipped, hydraulicdiamlist, rgaslist = rlist, fincoolingfactorfunc=fincoolingfactorfunc, dxlist = dxlist)

    material = "inconel 715"
    Structure = CS.StructuralAnalysis(rlist, xlist, nlist, chlist, cwlist, twlist, material)
    FOSlist = Structure.FOS(Twglist,Twclist,coolantpressurelist,preslist)
    xlist=np.flip(xlist) #idk whats flipping this haha but its something in the steadystatetemps function, so we have to flip it back
    return alistflipped, xlist, vlistflipped,\
     hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist,\
        Twglist, hglist, Qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist,\
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
        'throat_radius_curvature': .0254 *2,
        'dp': 150 * const.psiToPa,
        'impulseguess' :  495555.24828424345,
        'rc' : .11,
        'thetac' : (35*math.pi/180)} #623919}
configtitle = "NitroEthanol75 at 300 psi WorstCaseMass"
params = Main.main(args, configtitle, output=False)

params['thetac'] = (35*math.pi/180)
#conevol = math.pi*params['rc']**3*math.tan(params['thetac'])/3 - math.pi*params['rt']**3*math.tan(params['thetac'])/3
volfunc = lambda lc : math.pi*params['rc']**2*lc  +\
     math.pi*params['rc']**3/math.tan(params['thetac'])/3 -\
         math.pi*params['rt']**3/math.tan(params['thetac'])/3
lstarminimizer = lambda lc : volfunc(lc)/(params['rt']**2*math.pi) - params['lstar']
result = scipy.optimize.root(lstarminimizer, .05, args=(), method='hybr', jac=None, tol=None, callback=None, options=None)
params['lc']=result['x'][0]
xlist = np.linspace(0, params['lc'] + (params['rc'] - params['rt']) / math.tan(params['thetac']) + params['ln_conical'], 120)
       
rlist,xlist = RListGenerator.paraRlist(xlist, params['lc'], params['rc'],
                                  params['lc'] + (params['rc'] - params['rt'])/(math.tan(params['thetac'])),
                                  params['rt'],
                                  params['lc'] + (params['rc'] - params['rt'])/(math.tan(params['thetac'])) + params['ln_conical'],
                                  params['re'], params['lc']*1.8, 1.5*2*.0254, .4*2*.0254, math.pi/6 , 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re
# xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar

TC = ThrustChamber.ThrustChamber(rlist, xlist)
print(TC.rt, TC.xt, TC.xns)

machlist, preslist, templist = TC.flowSimple(params)


xlistflipped = np.flip(xlist)
rlistflipped  = np.flip(rlist)
chlist  = np.ones(len(xlist))*.003 #(TC.rt/rlistflipped)**.5
twlist  = (rlistflipped/TC.rt)*.001 
nlist  = np.ones(len(xlist))*80
ewlist  = chlist*2
#HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
helicitylist  = (rlistflipped**math.log(2,params['rc']/params['rt'])/TC.rt**math.log(2,params['rc']/params['rt']))*45*math.pi/180
chanelToLandRatio = 2
for index in range(0,np.size(helicitylist )):
    if helicitylist [index]>math.pi/2:
        helicitylist [index] = math.pi/2
for index in range(0,np.size(ewlist )):
    if ewlist [index]>.008:
        ewlist [index] = .008

alistflipped, xlist, vlistflipped,\
        hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist,\
            Twglist, hglist, Qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist, FOSlist = RunCoolingSystem(chlist,
    twlist,
    nlist,
    helicitylist,
    params,
    xlist,
    rlist,
    chanelToLandRatio, 
    TC )

plt.figure()

plt.plot(xlistflipped,FOSlist)
plt.xlabel("TC position, meters")
plt.ylabel("Factor of Safety")
plt.legend()
plt.title("Factor of Safety")

plt.figure()

plt.plot(xlistflipped,Twglist, label = "Gas Side Wall Temp [K]")
plt.plot(xlistflipped,Twclist, label = "Coolant Side Wall Temp [K]")
plt.plot(xlistflipped,Tclist, label = "Coolant Temp [K]")
plt.xlabel("TC position, meters")
plt.ylabel("Temp [K]")
plt.legend()
dxtrunc=(xlistflipped[0]-xlistflipped[1])
#plt.title(f"Gas Side Wall Temps, steps ={len(xlistflippedmat[0,:])}, dx = {'%.5f'%dxtrunc}")
plt.title(f"Temperatures")

plt.figure()

plt.plot(xlistflipped,Twclist)
plt.xlabel("TC position, meters")
plt.ylabel("Coolant SideWall Gas Temp [K]")
plt.legend()
plt.title("Coolant Side Wall Temp")

plt.figure()

plt.plot(xlistflipped,fincoolingfactorlist)
plt.xlabel("TC position, meters")
plt.ylabel('Fin Factor (Corrected coolant side area/gas side area)')
plt.legend()
plt.title("Fin Factor")

plt.figure()

plt.plot(xlistflipped,Tclist)
plt.xlabel("TC position, meters")
plt.ylabel("Coolant Temp [K]")
plt.legend()
plt.title("Coolant Temp")

plt.figure()

plt.plot(xlistflipped,coolantpressurelist/const.psiToPa)
plt.xlabel("TC position, meters")
plt.ylabel("Coolant pressure [psi]")
plt.legend()
plt.title("Coolant pressure")

plt.figure()

plt.plot(np.hstack((xlistflipped,xlist)),np.hstack((rlistflipped,-rlist)),'k')
plt.xlabel("TC position [m]")
plt.ylabel("Radius [m]")
plt.title("Thrust Chamber Shape")



plt.figure()
plt.plot(xlistflipped,chlist*1000,'r',label="Chanel Height [mm]")
plt.plot(xlistflipped,cwlist*1000,'b',label="Chanel Width [mm]")
plt.plot(xlistflipped,twlist*1000,'k',label="Wall Thickness [mm]")
plt.plot(xlistflipped,hydraulicdiamlist*1000,'g',label="Hydraulic Diam [mm]")
plt.plot(xlistflipped,helicitylist*180/math.pi/10,'m',label="helicity [10's of degrees]")
plt.plot(xlistflipped,vlistflipped,'c',label="Coolant Velocity [m/s]")
plt.legend()
plt.xlabel("TC position, meters")
plt.ylabel("Thicknesses, [mm]")
plt.title(f"Geometry , n = {int(nlist[1])}, Chanel to Land = {chanelToLandRatio}")

#plt.show()


pltxlist=np.flip(xlist) #idk whats flipping this haha but its something in the steadystatetemps function, so we have to flip it back
#xlist = np.flip(xlist) # flip it back damnit!
helicitylist = np.flip(helicitylist)

xlistnew, ylistnew,zlistnew,hydrualicdiamlistnew,chlistnew = CAD.ChanelBean(xlist,rlist,np.flip(twlist),
                helicitylist, np.flip(chlist), np.flip(cwlist)) 
plt.figure()

plt.plot(xlistnew[:,0],zlistnew[:,0])
plt.title("zig-zag bean lol")
path=os.path.join( "Configs","CAD_CDR_Engine_Bean_refined")
os.makedirs(path,exist_ok=True)
np.savetxt(os.path.join( path,"rlistflipped.csv"), rlistflipped, delimiter=",")   
np.savetxt(os.path.join( path,"twlist.csv"), twlist, delimiter=",") 
np.savetxt(os.path.join( path,"chlist.csv"), chlist, delimiter=",") 
np.savetxt(os.path.join( path,"cwlist.csv"), cwlist, delimiter=",") 
np.savetxt(os.path.join( path,"hydraulicdiamlist.csv"), hydraulicdiamlist, delimiter=",") 
np.savetxt(os.path.join( path,"nlist.csv"), nlist, delimiter=",") 
#with open(os.path.join(path, "chanelsweepcurve.sldcrv"), "w") as f:
#    for i in range(len(xchanel)):
#        #print(f"{xchanel[i]} {ychanel[i]} {zchanel[i]}", file=f) # this if for olivers axis's
#        print(f"{ychanel[i]} {xchanel[i]} {zchanel[i]}", file=f) # this if for roberts axis's
#with open(os.path.join(path, "chanelguidingcurve_height.sldcrv"), "w") as f:
#    for i in range(len(xchanel)):
#        #print(f"{xchanel[i]} {ychanel[i]} {zchanel[i]}", file=f) # this if for olivers axis's
#        print(f"{yheight[i]} {xheight[i]} {zheight[i]}", file=f) # this if for roberts axis's
#with open(os.path.join(path, "chanelguidingcurve_width.sldcrv"), "w") as f:
#    for i in range(len(xchanel)):
#        #print(f"{xchanel[i]} {ychanel[i]} {zchanel[i]}", file=f) # this if for olivers axis's
#        print(f"{ywidth[i]} {xwidth[i]} {zwidth[i]}", file=f) # this if for roberts axis's
for index in range(0,xlistnew.shape[0]):
    with open(os.path.join(path, f"chanelcurve{index}.sldcrv"), "w") as f:
        for i in range(len(xlistnew[1,:])):
            #print(f"{xchanel[i]} {ychanel[i]} {zchanel[i]}", file=f) # this if for olivers axis's
            print(f"{ylistnew[index,i]} {xlistnew[index,i]} {zlistnew[index,i]}", file=f) # this if for roberts axis's

with open(os.path.join(path, "internalradius.sldcrv"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {rlist[i]} {0}", file=f)
newxlist, externalrlist = CAD.rlistExtender(xlist,rlist,ewlist+np.flip(chlist)+np.flip(twlist))
with open(os.path.join(path, "externalradius.sldcrv"), "w") as f:
    for i in range(len(xlist)):
        print(f"{newxlist[i]} {externalrlist[i]} {0}", file=f)
with open(os.path.join(path, "midprofile.sldcrv"), "w") as f:
    for i in range(int(xlistnew.shape[0]/2)):
        print(f"{ylistnew[i*2,int(len(xlist)/2)]} {xlistnew[i*2,int(len(xlist)/2)]} {zlistnew[i*2,int(len(xlist)/2)]}", file=f) # this if for roberts axis's
    for i in np.flip(np.arange(int(xlistnew.shape[0]/2))):
        print(f"{ylistnew[i*2+1,int(len(xlist)/2)]} {xlistnew[i*2+1,int(len(xlist)/2)]} {zlistnew[i*2+1,int(len(xlist)/2)]}", file=f)
with open(os.path.join(path, "quarterprofile.sldcrv"), "w") as f:
    for i in range(int(xlistnew.shape[0]/2)):
        print(f"{ylistnew[i*2,int(len(xlist)/4)]} {xlistnew[i*2,int(len(xlist)/4)]} {zlistnew[i*2,int(len(xlist)/4)]}", file=f) # this if for roberts axis's
    for i in np.flip(np.arange(int(xlistnew.shape[0]/2))):
        print(f"{ylistnew[i*2+1,int(len(xlist)/4)]} {xlistnew[i*2+1,int(len(xlist)/4)]} {zlistnew[i*2+1,int(len(xlist)/4)]}", file=f)
with open(os.path.join(path, "threeqprofile.sldcrv"), "w") as f:
    for i in range(int(xlistnew.shape[0]/2)):
        print(f"{ylistnew[i*2,int(len(xlist)/4*3)]} {xlistnew[i*2,int(len(xlist)/4*3)]} {zlistnew[i*2,int(len(xlist)/4*3)]}", file=f) # this if for roberts axis's
    for i in np.flip(np.arange(int(xlistnew.shape[0]/2))):
        print(f"{ylistnew[i*2+1,int(len(xlist)/4*3)]} {xlistnew[i*2+1,int(len(xlist)/4*3)]} {zlistnew[i*2+1,int(len(xlist)/4*3)]}", file=f)

#NOW MAKE BOX CURVES

xlistnew, ylistnew,zlistnew = CAD.ChanelBoxCorners(xlist,rlist,np.flip(twlist),
                helicitylist, np.flip(chlist), np.flip(cwlist)) 
path=os.path.join( "Configs","CAD_CDR_Engine_Box_flippedtw_refined")
os.makedirs(path,exist_ok=True)

for index in range(0,4):
    with open(os.path.join(path, f"chanelcurve{index}.sldcrv"), "w") as f:
        for i in range(len(xlistnew[1,:])):
            #print(f"{xchanel[i]} {ychanel[i]} {zchanel[i]}", file=f) # this if for olivers axis's
            print(f"{ylistnew[index,i]} {xlistnew[index,i]} {zlistnew[index,i]}", file=f) # this if for roberts axis's

with open(os.path.join(path, "internalradius.sldcrv"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {rlist[i]} {0}", file=f)
newxlist, externalrlist = CAD.rlistExtender(xlist,rlist,ewlist+np.flip(chlist)+np.flip(twlist))
with open(os.path.join(path, "externalradius.sldcrv"), "w") as f:
    for i in range(len(xlist)):
        print(f"{newxlist[i]} {externalrlist[i]} {0}", file=f)
with open(os.path.join(path, "Twglist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {Twglist[i]}", file=f)
with open(os.path.join(path, "Twclist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {Twclist[i]}", file=f)
with open(os.path.join(path, "QdotTotallist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {Qdotlist[i]*nlist[i]*(abs(xlist[0]-xlist[1]))}", file=f)
with open(os.path.join(path, "qdotlist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {Qdotlist[i]*nlist[i]/(math.pi*rlist[i]*2)}", file=f)
with open(os.path.join(path, "QdotPerChanellist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {Qdotlist[i]*(abs(xlist[0]-xlist[1]))}", file=f)
with open(os.path.join(path, "Machlist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {machlist[i]}", file=f)
with open(os.path.join(path, "GasTemplist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {templist[i]}", file=f)
with open(os.path.join(path, "GasPressurelist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {preslist[i]}", file=f)
with open(os.path.join(path, "CoolantPressurelist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {coolantpressurelist[i]}", file=f)
with open(os.path.join(path, "CoolantTemperaturelist.txt"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {Tclist[i]}", file=f)


thetalist = np.arange(0, 2*math.pi, .1)
theta, r = np.meshgrid(thetalist, rlist-.01)
theta, xgrid = np.meshgrid(thetalist, xlist)
zgrid = r*np.cos(thetalist)
ygrid = r*np.sin(thetalist)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

for index in range(0,xlistnew.shape[0]):
    ax.plot(xlistnew[index,], zlistnew[index,], ylistnew[index,],linewidth=.5)
#chanel0 = ax.plot(xlistnew[0,], zlistnew[0,], ylistnew[0,],'r',linewidth=.5)
#chanel1 = ax.plot(xlistnew[1,], zlistnew[1,], ylistnew[1,],'g',linewidth=.5)
#chanel2 = ax.plot(xlistnew[2,], zlistnew[2,], ylistnew[2,],'b',linewidth=.5)
#chanel3 = ax.plot(xlistnew[3,], zlistnew[3,], ylistnew[3,],'m',linewidth=.5)
surf = ax.plot_surface(zgrid,ygrid,xgrid, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True, alpha=.5)

# Customize the z axis.
ax.set_zlim(0, np.max(xlist))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)


plt.show()

print("end")