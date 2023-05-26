import sys
sys.path.insert(1,"./")
import numpy as np
import sympy as sym
import math
import difflib
from rocketcea.cea_obj_w_units import CEA_Obj
import Toolbox.RocketCEAAssister as RA
import Toolbox.IsentropicEquations as IE
from rocketprops.rocket_prop import get_prop
import re as regex
import Toolbox.Constant as const
from scipy import interpolate
import matplotlib.pyplot as plt
def optimalMr(args, plot = False):
    try:
        args['pc']
    except:
        print("YOU SHOULD REALLY INPUT A PC, JUST ASSUMING 700 psi NOW")
        args['pc']=700*const.psiToPa
    pambient=14.7*const.psiToPa
    pc=args['pc']


    try:
        RA.makeEthanolBlend(int(regex.search(r'\d+', args['fuelname']).group()))
    except:
        print(args['fuelname'])
    try: #maybe there is a kyename error
        if args['cr'] is None:
            CEA = CEA_Obj(oxName=args['oxname'], fuelName=args['fuelname'], useFastLookup=1, isp_units='sec',
                                    cstar_units='m/s',
                                    pressure_units='Pa', temperature_units='K', sonic_velocity_units='m/s',
                                    enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units='J/kg-K',
                                    viscosity_units='millipoise', thermal_cond_units='W/cm-degC', fac_CR=None,
                                    make_debug_prints=False)
        else:
            CEA = CEA_Obj(oxName=args['oxname'], fuelName=args['fuelname'], useFastLookup=1, isp_units='sec',
                                    cstar_units='m/s',
                                    pressure_units='Pa', temperature_units='K', sonic_velocity_units='m/s',
                                    enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units='J/kg-K',
                                    viscosity_units='millipoise', thermal_cond_units='W/cm-degC', fac_CR=args['cr'],
                                    make_debug_prints=False)
    except KeyError:
        CEA = CEA_Obj(oxName=args['oxname'], fuelName=args['fuelname'], useFastLookup=1, isp_units='sec',
                                cstar_units='m/s',
                                pressure_units='Pa', temperature_units='K', sonic_velocity_units='m/s',
                                enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units='J/kg-K',
                                viscosity_units='millipoise', thermal_cond_units='W/cm-degC', fac_CR=None,
                                make_debug_prints=False)
    philist = np.arange(.5,4,.01)
    mrlist = np.zeros( philist.size)
    isplistfrozen = np.zeros( philist.size)
    isplist = np.zeros( philist.size)
    epslistfrozen = np.zeros( philist.size)
    epslist = np.zeros( philist.size)
    ind=0
    for phi in philist:
        mrlist[ind]=CEA.getMRforER(ERphi=philist[ind])
        epslist[ind] = CEA.get_eps_at_PcOvPe(Pc=pc, MR=mrlist[ind], PcOvPe=pc/pambient, frozen=0, frozenAtThroat=0)
        epslistfrozen[ind] = CEA.get_eps_at_PcOvPe(Pc=pc, MR=mrlist[ind], PcOvPe=pc / pambient, frozen=1,
                                                       frozenAtThroat=1)
        isplist[ind]= CEA.estimate_Ambient_Isp(Pc=pc, MR=mrlist[ind], eps=epslist[ind],
                                                           Pamb=pambient,
                                                           frozen=0, frozenAtThroat=0)[0]
        isplistfrozen[ind] = CEA.estimate_Ambient_Isp(Pc=pc, MR=mrlist[ind], eps=epslist[ind],
                                                           Pamb=pambient,
                                                           frozen=1, frozenAtThroat=1)[0]
        ind = ind+1



    ispmaxeq = isplist.max()
    ispmaxfrozen = isplistfrozen.max()
    ispavglist = np.array((isplist + isplistfrozen) / 2)
    ispmaxavg = ispavglist.max()
    mrideal = mrlist[ispavglist.argmax()]
    phiideal = philist[ispavglist.argmax()]
    if plot:
        title = f"Optimal MR Plots for {args['fuelname']} and {args['oxname']} at Pc = {args['pc'] / const.psiToPa}"
        fig, axs = plt.subplots(3, 1)
        fig.suptitle(title)

        axs[0].plot(philist, isplist, 'g', label="isp equilibrium")  # row=0, column=0
        axs[0].plot(philist, isplistfrozen, 'r', label="isp frozen")  # row=0, column=0
        axs[0].plot(philist, (isplist+isplistfrozen)/2, 'm', label="isp average")  # row=0, column=0
        axs[0].plot([philist[isplist.argmax()],philist[isplist.argmax()]],[0, 1000], 'k--')
        axs[0].plot( [philist[isplistfrozen.argmax()], philist[isplistfrozen.argmax()]],[0, 1000], 'k--')
        axs[0].plot([philist[ispavglist.argmax()], philist[ispavglist.argmax()]],[0, 1000], 'm--')
        axs[1].plot(mrlist, isplist, 'g', label="isp equilibrium")  # row=0, column=0
        axs[1].plot(mrlist, isplistfrozen, 'r', label="isp frozen")  # row=0, column=0
        axs[1].plot(mrlist, (isplist + isplistfrozen) / 2, 'm', label="isp average")  # row=0, column=0
        axs[1].plot( [mrlist[isplist.argmax()], mrlist[isplist.argmax()]],[0, 1000], 'k--')
        axs[1].plot( [mrlist[isplistfrozen.argmax()], mrlist[isplistfrozen.argmax()]],[0, 1000], 'k--')
        axs[1].plot( [mrlist[ispavglist.argmax()], mrlist[ispavglist.argmax()]],[0, 1000], 'm--')
        axs[2].plot(mrlist, (epslist + epslistfrozen) / 2, 'g', label="isp average")

        axs[0].set_title(f"Isp Vs Phi ")
        axs[1].set_title(f"Isp Vs Mr ")
        axs[2].set_title(f"Supersonic Area Ratio Vs Mr")
        axs[0].legend()
        axs[1].legend()
        axs[0].set(ylabel='Isp')
        axs[1].set(ylabel='Isp')
        axs[1].set_xlim([0,mrlist.max()])
        axs[0].set_xlim([0, philist.max()])
        axs[1].set_ylim([150, ispavglist.max()+30])
        axs[0].set_ylim([150, ispavglist.max()+30])

    return ispmaxavg, mrideal, phiideal, ispmaxeq, ispmaxfrozen
