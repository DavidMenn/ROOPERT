#from Components.ThrustChamber import ThrustChamber
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant
import numpy as np
import Toolbox.RListGenerator as RListGenerator
import Toolbox.RocketCEAAssister as RA
import matplotlib.pyplot as plt

xlist = np.linspace(0, .35, 1000)
rlist = RListGenerator.sharpRList(xlist, .2, .07, .24, .02, .35, .06) #xlist, xns, rc, xt, rt, xe, re
#TC = ThrustChamber(rlist,xlist)
#print(TC.rt, TC.xt, TC.xns)

RA.makeEthanolBlend(75)

C = CEA_Obj(oxName="LOX", fuelName="Ethanol_75")

s = C.get_full_cea_output( Pc=350, MR=[1.19,1.2,1.21] ,eps=4, short_output=1)

print( s )

