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
import Main
from scipy.optimize import minimize_scalar
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




#find rli
#find rlist
# preset geometry
# run temps
# run structures
# return FS's



"""args = {
    'thrust': 17010.52372
,  # Newtons
    'time': 31.48,  # s
    # 'rho_ox' : 1141, #Kg/M^3
    # 'rho_fuel' : 842,
    'pc': 300 * const.psiToPa,
    'pe': 14.7 * const.psiToPa,
    'phi': 1.089,
    # TYPICALLY WE SET 'phi'
    'cr': 4,
    'lstar': 1,
    'fuelname': 'Ethanol',
    'oxname': 'LOX',
    'throat_radius_curvature': .0254 / 2,
    'dp': 150 * const.psiToPa}"""
#params = FAC.SpreadsheetSolver(args)

args = {
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 10 * const.psiToPa,
       # 'phi':1,
        'cr': 4,
        'TWR' : 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa,
        'impulseguess' : 494381.0297091524} #623919}
configtitle = "NitroEthanol75 at 300 psi Optimize_Wall_Thickness"
params = Main.main(args, configtitle, output=False)


#rlist = RListGenerator.sharpRList(xlist, params['lc'], params['rc'],
#                                  params['lc'] + (params['rc'] - params['rt']) / (2 ** .5),
#                                  params['rt'],
#                                  params['lc'] + (params['rc'] - params['rt']) / (2 ** .5) + params['ln_conical'],
#                                  params['re'])  #xlist, xns, rc, xt, rt, xe, re, rcf, rtf, theta, phi

#rlist,xlist = RListGenerator.roundRList(xlist, params['lc'], params['rc'],
#                                  params['lc'] + (params['rc'] - params['rt']) / (2 ** .5),
#                                  params['rt'],
#                                  params['lc'] + (params['rc'] - params['rt']) / (2 ** .5) + params['ln_conical'],
#                                  params['re'], .054, .0254)  # xlist, xns, rc, xt, rt, xe, re

#plt.plot(xlist,rlist)
xlist = np.linspace(0, params['lc'] + (params['rc'] - params['rt']) / (2 ** .5) + params['ln_conical'], 500)
params['lc']=params['lc']*.5
#paraRlist (xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar)
rlist,xlist = RListGenerator.paraRlist(xlist, params['lc'], params['rc'],
                                  params['lc'] + (params['rc'] - params['rt'])/(math.tan(20*math.pi/180)),
                                  params['rt'],
                                  params['lc'] + (params['rc'] - params['rt'])/(math.tan(20*math.pi/180)) + params['ln_conical'],
                                  params['re'], params['lc']*1.8, 3*.0254, 3*.0254, math.pi/6, 8*math.pi/180, params['er'])  # xlist, xns, rc, xt, rt, xe, re

#plt.plot(paraxlist,pararlist)
#plt.show()
TC = ThrustChamber.ThrustChamber(rlist, xlist)
print(TC.rt, TC.xt, TC.xns)

machlist, preslist, templist = TC.flowSimple(params)


xlistflipped = np.flip(xlist)
rlistflipped = np.flip(rlist)
chlist = (TC.rt/rlistflipped)**.5*.003 
twlist = (rlistflipped/TC.rt)*.001 
nlist = np.ones(len(xlist))*40
ewlist = np.ones(len(xlist))*.005
#HELICITY IS DEFINED AS 90 DEGREES BEING A STAIGHT CHANEL, 0 DEGREES BEING COMPLETILY CIRCUMFRNEITAL
helicitylist = (rlistflipped**1.5/TC.rt**1.5)*30*math.pi/180

for index in range(0,np.size(helicitylist)):
    if helicitylist[index]>math.pi/2:
        helicitylist[index] = math.pi/2

alistflipped, nlist, coolingfactorlist, heatingfactorlist, xlist, vlistflipped, twlistflipped,\
     hydraulicdiamlist, salistflipped, dxlist, fincoolingfactorfunc, cwlist   = fixedRSquareChanelSetup(params = params,
                                        xlist = xlistflipped, rlist = rlistflipped,chlist = chlist ,
                                        chanelToLandRatio = 2 ,twlist = twlist ,nlist = nlist,
                                        helicitylist=helicitylist)



 # HAVE THIS COMMETNED OUT WHILE IM MAKING GEOMS TO AVOID HAVING TO CALC TEMPS EVERY TIME

Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, fincoolingfactorlist, rholist, viscositylist, Relist = CS.steadyStateTemperatures(
    None, TC, params, salistflipped, nlist, coolingfactorlist,
    heatingfactorlist, xlistflipped, vlistflipped, 293, params['pc'] + params['pc']*.2 + 50*const.psiToPa, twlistflipped, hydraulicdiamlist, rgaslist = rlist, fincoolingfactorfunc=fincoolingfactorfunc, dxlist = dxlist)
#def steadyStateTemperatures(wallmaterial, thrustchamber, params, salist, nlist, coolingfactorlist,
 #                           heatingfactorlist, xlist, vlist, initialcoolanttemp, initialcoolantpressure, twlist,
 #                           hydraulicdiamlist, rgaslist = None, fincoolingfactorfunc = None, dxlist = None, helicitylist = None):
material = "inconel 715"
Structure = CS.StructuralAnalysis(rlist, xlist, nlist, chlist, cwlist, twlist, material)
FOSlist = Structure.FOS(Twglist,Twclist,coolantpressurelist,preslist)

xlist=np.flip(xlist) #idk whats flipping this haha but its something in the steadystatetemps function, so we have to flip it back

helicitylist = np.flip(helicitylist)

#xlistnew, ylistnew,zlistnew = CAD.ChanelBoxCorners(xlist,rlist,twlist,
#                helicitylist, np.flip(chlist), np.flip(cwlist)) 
xlistnew, ylistnew,zlistnew = CAD.ChanelBean(xlist,rlist,twlist,
                helicitylist, np.flip(chlist), np.flip(cwlist)) 

path=os.path.join( "Configs","CAD_Curve")
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
        for i in range(len(xlist)):
            #print(f"{xchanel[i]} {ychanel[i]} {zchanel[i]}", file=f) # this if for olivers axis's
            print(f"{ylistnew[index,i]} {xlistnew[index,i]} {zlistnew[index,i]}", file=f) # this if for roberts axis's
with open(os.path.join(path, "internalradius.sldcrv"), "w") as f:
    for i in range(len(xlist)):
        print(f"{xlist[i]} {rlist[i]} {0}", file=f)
newxlist, externalrlist = CAD.rlistExtender(xlist,rlist,ewlist+np.flip(chlist)+np.flip(twlist))
with open(os.path.join(path, "externalradius.sldcrv"), "w") as f:
    for i in range(len(xlist)):
        print(f"{newxlist[i]} {externalrlist[i]} {0}", file=f)


thetalist = np.arange(0, 2*math.pi, .1)
theta, r = np.meshgrid(thetalist, rlist-.01)
theta, xgrid = np.meshgrid(thetalist, xlist)
zgrid = r*np.cos(thetalist)
ygrid = r*np.sin(thetalist)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

for index in range(0,xlistnew.shape[0]):
    ax.plot(xlistnew[index,], zlistnew[index,], ylistnew[index,],linewidth=.5)
surf = ax.plot_surface(zgrid,ygrid,xgrid, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True, alpha=.5)

# Customize the z axis.
ax.set_zlim(0, xlist[-1])
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.show()

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
plt.plot(xlistflipped,chlist,'r',label="Chanel Height(mm)")
plt.plot(xlistflipped,cwlist,'b',label="Chanel Widht(mm)")
plt.plot(xlistflipped,twlist,'k',label="Wall Thickness(mm)")
plt.plot(xlistflipped,hydraulicdiamlist,'g',label="Hydraulic Diam(mm)")
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
