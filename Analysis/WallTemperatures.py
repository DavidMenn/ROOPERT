import sys
sys.path.insert(1,"./")
import scipy.optimize
#from Components.ThrustChamber import ThrustChamber
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer, add_new_propellant
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
import Analysis.FirstOrderCalcs as FAC
import Components.ThrustChamber as ThrustChamber
import Components.CoolingSystem as CS
import Components.StructuralApproximation as SA
from scipy.optimize import minimize_scalar

"""THIS IS SUPER AIDS BUT FOR SOME REASON TANSPORt PROPS ARE KINDA BORKEN IN ROCKETCEA OUT OF NOWHERE
I change fast look up to 1 and that somehow fixed it but this could be some long term problem!
here is code to get the full output and just compare:

from rocketcea import cea_obj
C=cea_obj.CEA_Obj( oxName='LOX', fuelName='Ethanol', fac_CR=4,
        useFastLookup=0, # deprecated
        makeOutput=0, make_debug_prints=True) 
out=C.get_full_cea_output(Pc=100.0, MR=1.0, eps=40.0, subar=4, PcOvPe=100/14.7,
                            frozen=0, frozenAtThroat=0, short_output=0, show_transport=1,
                            pc_units='psia', output='calories', show_mass_frac=False,
                            fac_CR=4)
                           
                            """
args = {
        'thrust': 1200 * const.lbToN,  # Newtons
        'time': 7.5,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 700 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi' : 1,
    #TYPICALLY WE SET 'phi'
        'cr' :4,
        'lstar' : 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254/2,
        'dp' : 100*const.psiToPa}
args = {
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 700 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi' : 1,
    #TYPICALLY WE SET 'phi'
        'cr' :4,
        'lstar' : 1,
        'fuelname': 'Ethanol',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254/2,
        'dp' : 100*const.psiToPa}
params = FAC.SpreadsheetSolver(args)
print(f"of ratio is {params['rm']}")
xlist = np.linspace(0, params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], 1000)

rlist = RListGenerator.sharpRList(xlist, params['lc'], params['rc'], params['lc']+(params['rc']-params['rt'])/(2**.5),
                                  params['rt'], params['lc']+(params['rc']-params['rt'])/(2**.5)+params['ln_conical'], params['re']) #xlist, xns, rc, xt, rt, xe, re
TC = ThrustChamber.ThrustChamber(rlist,xlist)
print(TC.rt, TC.xt, TC.xns)

machlist,preslist,templist = TC.flowSimple(params)







chanelthickness = (rlist==TC.rc)*.03*.0254+(rlist!=TC.rc)*.03*.0254
tw=.1*.0254
alist=math.pi*((rlist+tw+chanelthickness)**2-(rlist+tw)**2)
alistflipped, n, coolingfactorlist, heatingfactorlist, xlistflipped, vlistflipped, twlistflipped, hydraulicdiamlist, salistflipped = CS.coaxialShellSetup(TC,
                                        params,rlist, tw=tw, chanelthickness = chanelthickness, helicity=1, dt=.002,vlist = None, alist=alist)

Twglist, hglist, qdotlist, Twclist, hclist, Tclist, coolantpressurelist, qdotlist, Trlist, rholist, viscositylist, Relist = CS.steadyStateTemperatures(None,TC, params, salistflipped,n, coolingfactorlist,
                            heatingfactorlist, xlistflipped, vlistflipped ,293, params['pc']+params['dp'], twlistflipped, hydraulicdiamlist)


title="testing flowSimple for thrustchamber"
fig, axs = plt.subplots(4,3)
fig.suptitle(title)

axs[0,1].plot(xlistflipped,hglist , 'g')  # row=0, column=0
axs[1,1].plot(xlistflipped,hclist , 'r')# row=1, column=0
axs[2,1].plot(xlistflipped,Trlist, 'b') # row=0, column=0
#axs[3,1].plot(np.hstack((np.flip(TC.xlist),xlist)),np.hstack((np.flip(rlist),-rlist)) , 'k') # row=0, column=0
axs[3,1].plot(xlistflipped,qdotlist, 'r')# row=1, column=0

axs[0,1].set_title('hglist')
axs[1,1].set_title('hclist')
axs[2,1].set_title('trlist')
axs[3,1].set_title('qdot')

axs[0,0].plot(xlistflipped,Twglist , 'g', label="Gas Side Wall Temp")
axs[0,0].plot(xlistflipped,Twclist , 'r', label="CoolantSide Wall Temp") # row=0, column=0
axs[0,0].plot(xlistflipped,Tclist , 'b', label="Coolant Temp") #
axs[1,0].plot(xlistflipped,Tclist , 'r')# row=1, column=0
axs[2,0].plot(xlistflipped,hydraulicdiamlist , 'r')# row=1, column=0
axs[3,0].plot(xlistflipped,coolantpressurelist /const.psiToPa, 'k') # row=0, column=0

axs[0,0].set_title('Twg')
axs[1,0].set_title('Tc')
axs[2,0].set_title('hydraulicdiam')
axs[3,0].set_title('coolantpressure (psi)')
axs[0,0].legend()

axs[0,2].plot(xlistflipped,Twglist*const.degKtoR-458.67 , 'g', label="Gas Side Wall Temp, F")
#axs[0,2].plot(xlistflipped,Tcoatinglist*const.degKtoR-458.67 , 'k', label="Opposite side of coating Temp, F")
axs[0,2].plot(xlistflipped,Twclist * const.degKtoR-458.67, 'r', label="CoolantSide Wall Temp, F") # row=0, column=0

axs[0,2].plot(xlistflipped,Tclist * const.degKtoR-458.67, 'b', label="Coolant Temp, F") #
axs[1,2].plot(xlistflipped,Relist, 'r')# row=1, column=0
axs[2,2].plot(xlistflipped,salistflipped , 'r')# row=1, column=0
axs[3,2].plot(xlistflipped,rholist, 'k') # row=0, column=0

axs[0,2].set_title('Twg')
axs[1,2].set_title('Reynolds Number of Coolant')
axs[2,2].set_title('salist')
axs[3,2].set_title('density of coolant')
axs[0,2].legend()

print(f"max twg = {np.max(Twglist)} in kelvin, {np.max(Twglist)*const.degKtoR} in Rankine (freedom)\n max Twc ="
      f" {np.max(Twclist)} in kelvin, {np.max(Twclist)*const.degKtoR} in Rankine (freedom)")
# Hide x labels and tick labels for top plots and y ticks for right plots.


plt.show()


import Toolbox.Constant as const
import Analysis.FirstOrderCalcs as FAC
args = {
'pc': 700 * const.psiToPa,
'phi' : 1,
'fuelname': 'Ethanol',
'oxname': 'N2O'}
params = FAC.SpreadsheetSolver(args)