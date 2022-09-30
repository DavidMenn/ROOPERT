"""David Menn
The intent of this class is to model everything that goes on inside the combustion chamber and nozzle.
A ThrustChamber object just stores the geometry of the combustion chamber, with xlist being the list
of x values (in m) for each array index corresponding to rlist, which is the radius of the crosssection
at that x value.

The utility of this class is that, once given input parameters from the injection, it can give you the
temp, pressure, axial velocity, and whatever else I decide of the mixture inside the combustion chamber,
which will then be stored in the object.

Should be able to spit out thrust, isp, mdot, etc.

Should have independent function that spits out data for the end of just CC portion to be used when
optimizing nozzle

Future additions will be using cantera to model the combustion and get a reasonable value for L*

Future additions will be the modelling of combustion instability and resonance (first order stuff)
"""
import numpy as np
#from rocketcea.cea_obj import CEA_Obj
import math
import Toolbox.IsentropicEquations as Ise

class ThrustChamber:

    def __init__(self, rlist, xlist):
        self.rlist = rlist
        self.xlist = xlist
        self.alist = math.pi*rlist**2

        self.rt = np.amin(rlist)
        self.xt = xlist[ np.argmin(rlist) ]
        self.at = math.pi * self.rt ** 2

        self.xns = xlist[ np.where( np.diff( rlist ) < 0)[0][0] ] #finds first instance of deviation, ns stands for nozzle start

        self.cr = (rlist[0]**2)/(self.rt**2)
        self.eps = (rlist[-1]**2)/(self.rt**2)




#Flow props are dicts in format {'name': String, 'Mdot' : int (kg/s),
""" Currently calculates like this: Assume known mdot, which is found by targetting an isp and thrust. Now we know from choked flow at the throat, assuming total temperature stays
constant throughout the nozzle, the total pressure at the throat and thus the pressure at the throat. Now we can simply work backwards to get the pressures and velocities
"""
    def flow(self, fuel_props, ox_props, pcEstimate = None): #MAKE SURE TO FIX UNITS OR YOURE GONNA GET SOME WEIRD SHIT

        convert from mdot/mdot to MR using getMRforER(ERphi=None, ERr=None)
        mdot = totalmdot
        CEA=CEA_Obj(propName='', oxName='', fuelName='', fac_CR=self.cr, useFastLookup=0, makeOutput=0, make_debug_prints=False) #maybe try to set units here
        if pcEstimate is None
            pcEstimate = 400 # idk some standard value i dont think this matters???

        initialTemp = get_Temperatures(Pc=pcEstimate, MR=mr, eps=40.0, frozen=0, frozenAtThroat=0)[1]
        initalProperties = CEA.get_Temperatures(Pc=pcEstimate, MR=mr, eps=40.0, frozen=0)
        initialPt = Ise.PFromMdotAtThroat(mdot , initalProperties[1] , Ise.totalT(initialTemp,initalProperties[1],1) , self.at , Ise.gasConstant()/initalProperties[0])
        get from that to Pc using area rule (put that in the isentopic toob ox)

        Do it again like three times

        Compare results