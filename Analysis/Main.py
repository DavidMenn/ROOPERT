from Components.ThrustChamber import ThrustChamber
import numpy as np
import Toolbox.RListGenerator as RListGenerator
import matplotlib.pyplot as plt
xlist = np.linspace(0, .35, 1000)
rlist = RListGenerator.sharpRList(xlist, .2, .07, .24, .02, .35, .06) #xlist, xns, rc, xt, rt, xe, re
TC = ThrustChamber(rlist,xlist)
print(TC.rt, TC.xt, TC.xns)
