import sys
sys.path.insert(1,"./")
import scipy.optimize
#from Components.ThrustChamber import ThrustChamber
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer, add_new_propellant
import numpy as np
import math
import Components.ThrustChamber as ThrustChamber
import Components.CoolingSystem as CS
from rocketcea.cea_obj_w_units import CEA_Obj
import Toolbox.RListGenerator as RListGenerator
import Toolbox.RocketCEAAssister as RA
import os
import Toolbox.IsentropicEquations as IE
import Toolbox.RocketEquation as RE
import difflib
import re as regex
from rocketprops.rocket_prop import get_prop
import Toolbox.Constant as const
import DetermineOptimalMR as DOMR
import matplotlib.pyplot as plt
import FirstOrderCalcs as FAC
import Components.ThrustChamber as ThrustChamber
import Components.StructuralApproximation as SA
from scipy.optimize import minimize_scalar
from Main import main
'''
args = {
        'thrust': 4000 * const.lbToN,  # Newtons, INITIAL GUESS
        'time': 40,  # s, INIT GUESS
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        # 'phi':1,
        'cr': 4,
        'TWR' : 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa}
configtitle="NitroEthanol75"
path=os.path.join( "Configs","Tradestudies",configtitle)
os.makedirs(path,exist_ok=False)


pressures = np.arange(150,750,25)*const.psiToPa
impulselist = np.zeros(pressures.size)
milist = np.zeros(pressures.size)
heightflist = np.zeros(pressures.size)
heightolist = np.zeros(pressures.size)
newheightlist=np.zeros(pressures.size)
thrustlist = np.zeros(pressures.size)
mdotlist = np.zeros(pressures.size)
mdotoxlist = np.zeros(pressures.size)
mdotfuellist = np.zeros(pressures.size)
isplist = np.zeros(pressures.size)
maxtwglist = np.zeros(pressures.size)
maxtwclist = np.zeros(pressures.size) #HAD TO COMMENT OUT LINE 317 and 318 incoling system!!
wstructlist = np.zeros(pressures.size)
lambdalist = np.zeros(pressures.size)
totalmasslist = np.zeros(pressures.size)
ind = 0
params=dict(impulse=500000)
for pres in pressures:
        args['pc']=pres
        print(pres)
        args['impulseguess']=params['impulse']
        params=main(args,"",output=False)
        impulselist[ind] = params['impulse']
        milist[ind] = params['mi']
        heightflist[ind] = params['heightfuel']
        heightolist[ind] = params['heightox']
        newheightlist[ind] = params['newheight']
        thrustlist[ind] = params['thrust']
        mdotlist[ind] = params['mdot']
        mdotoxlist[ind] = params['mdot_ox']
        mdotfuellist[ind] = params['mdot_fuel']
        isplist[ind] = params['isp']
        maxtwglist[ind] = params['twg_max']
        maxtwclist[ind] = params['twc_max']
        wstructlist[ind] = params['wstruct']
        lambdalist[ind] = params['L']
        totalmasslist[ind]=params['M']
        ind = ind+1


plt.plot(pressures/const.psiToPa,impulselist,'g')
plt.title("pressures (psi) vs impulse (ns)")
plt.savefig(os.path.join(path, "impulse.png"))
plt.figure()
plt.plot(pressures/const.psiToPa,isplist,'m')
plt.title("pressures (psi) vs isp (s)")
plt.savefig(os.path.join(path, "isp.png"))
plt.figure()
plt.plot(pressures/const.psiToPa,milist,'r')
plt.title("pressures (psi) vs inert mass (kg)")
plt.savefig(os.path.join(path, "inertmass.png"))
plt.figure()
plt.plot(pressures/const.psiToPa,thrustlist/const.lbToN)
plt.title("pressures (psi) vs thrust(lb)")
plt.savefig(os.path.join(path, "thrust.png"))
plt.figure()
plt.plot(pressures/const.psiToPa,newheightlist,'k')
plt.title("pressures (psi) vs heights (m)")
plt.savefig(os.path.join(path, "height.png"))
plt.figure()
plt.plot(pressures/const.psiToPa,mdotlist,'k', label='total')
plt.plot(pressures/const.psiToPa,mdotoxlist,'b', label='ox')
plt.plot(pressures/const.psiToPa,mdotfuellist,'r', label='fuel')
plt.title("pressures (psi) vs mass flow rate (kg/s)")
plt.legend()
plt.savefig(os.path.join(path, "mdot.png"))
plt.figure()

plt.plot(pressures/const.psiToPa,maxtwglist,'r', label='max twg')
plt.plot(pressures/const.psiToPa,maxtwclist,'b', label='max twc')
plt.title("pressures (psi) vs max wall temp (kelvin)")
plt.savefig(os.path.join(path, "chambertemps.png"))
plt.legend()

plt.show()

"""
#THIS IS THE TWR TRADE STUDY
args = {
        'thrust': 4000 * const.lbToN,  # Newtons, INITIAL GUESS
        'time': 40,  # s, INIT GUESS
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        # 'phi':1,
        'cr': 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'LOX',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa}
configtitle="LOXEthanol75_350psi_TWRstudy"
path=os.path.join("Configs","Tradestudies",configtitle)
os.makedirs(path,exist_ok=False)


#twrs = np.arange(150,751,25)*const.psiToPa
twrs = np.arange(1.25,5,.25) #THIS IS REALLY TWR
impulselist = np.zeros(twrs.size)
milist = np.zeros(twrs.size)
heightflist = np.zeros(twrs.size)
heightolist = np.zeros(twrs.size)
newheightlist=np.zeros(twrs.size)
thrustlist = np.zeros(twrs.size)
mdotlist = np.zeros(twrs.size)
mdotoxlist = np.zeros(twrs.size)
mdotfuellist = np.zeros(twrs.size)
isplist = np.zeros(twrs.size)
maxtwglist = np.zeros(twrs.size)
maxtwclist = np.zeros(twrs.size) #HAD TO COMMENT OUT LINE 317 and 318 incoling system!!
wstructlist = np.zeros(twrs.size)
lambdalist = np.zeros(twrs.size)
totalmasslist = np.zeros(twrs.size)
ind = 0
params=dict(impulse=500000)
for twr in twrs:
        args['TWR']=twr
        args['impulseguess']=params['impulse']
        print(twr)
        params=main(args,"",output=False)
        impulselist[ind] = params['impulse']
        milist[ind] = params['mi']
        heightflist[ind] = params['heightfuel']
        heightolist[ind] = params['heightox']
        newheightlist[ind] = params['newheight']
        thrustlist[ind] = params['thrust']
        mdotlist[ind] = params['mdot']
        mdotoxlist[ind] = params['mdot_ox']
        mdotfuellist[ind] = params['mdot_fuel']
        isplist[ind] = params['isp']
        maxtwglist[ind] = params['twg_max']
        maxtwclist[ind] = params['twc_max']
        wstructlist[ind] = params['wstruct']
        lambdalist[ind] = params['L']
        totalmasslist[ind]=params['M']
        ind = ind+1


plt.plot(twrs,impulselist,'g')
#plt.title("twrs (psi) vs impulse (ns)")
plt.title("TWR vs impulse (ns)")
plt.savefig(os.path.join(path, "impulse.png"))
plt.figure()
plt.plot(twrs,isplist,'m')
plt.title("TWR vs isp (s)")
plt.savefig(os.path.join(path, "isp.png"))
plt.figure()
plt.plot(twrs,milist,'r')
plt.title("TWR vs inert mass (kg)")
plt.savefig(os.path.join(path, "inertmass.png"))
plt.figure()
plt.plot(twrs,thrustlist/const.lbToN)
plt.title("TWR vs sea levelthrust(lb)")
plt.savefig(os.path.join(path, "thrust.png"))
plt.figure()
plt.plot(twrs,newheightlist,'k')
plt.title("TWR vs heights (m)")
plt.savefig(os.path.join(path, "height.png"))
plt.figure()
plt.plot(twrs,mdotlist,'k', label='total')
plt.plot(twrs,mdotoxlist,'b', label='ox')
plt.plot(twrs,mdotfuellist,'r', label='fuel')
plt.title("TWR vs mass flow rate (kg/s)")
plt.legend()
plt.savefig(os.path.join(path, "mdot.png"))
plt.figure()

plt.plot(twrs,maxtwglist,'r', label='max twg')
plt.plot(twrs,maxtwclist,'b', label='max twc')
plt.title("TWR vs max wall temp (kelvin)")
plt.savefig(os.path.join(path, "chambertemps.png"))
plt.legend()

plt.show()
"""
'''
"""
args = {
        'thrust': 4000 * const.lbToN,  # Newtons, INITIAL GUESS
        'time': 40,  # s, INIT GUESS
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        # 'phi':1,
        'TWR' : 4,
        'cr': 4,
        'lstar': 1,
        'fuelname': 'Ethanol_60',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa}
configtitle="NitroEthanol_300psi_waterpercentstudy_TWR4"
path=os.path.join("Configs","Tradestudies",configtitle)
os.makedirs(path,exist_ok=False)


#twrs = np.arange(150,751,25)*const.psiToPa
percents = np.array([60,65,70,75,80,85,90,95,99]) 
percenttitles= {60:'Ethanol_60',
65:'Ethanol_65',
70:'Ethanol_70',
75:'Ethanol_75',
80:'Ethanol_80',
85:'Ethanol_85',
90:'Ethanol_90',
95:'Ethanol_95',
99:'Ethanol_99',
}
impulselist = np.zeros(percents.size)
milist = np.zeros(percents.size)
heightflist = np.zeros(percents.size)
heightolist = np.zeros(percents.size)
newheightlist=np.zeros(percents.size)
thrustlist = np.zeros(percents.size)
mdotlist = np.zeros(percents.size)
mdotoxlist = np.zeros(percents.size)
mdotfuellist = np.zeros(percents.size)
isplist = np.zeros(percents.size)
maxtwglist = np.zeros(percents.size)
maxtwclist = np.zeros(percents.size) #HAD TO COMMENT OUT LINE 317 and 318 incoling system!!
wstructlist = np.zeros(percents.size)
lambdalist = np.zeros(percents.size)
totalmasslist = np.zeros(percents.size)
ind = 0
for perc in percents:
        args['fuelname']=percenttitles[perc]
        print(percenttitles[perc])
        params=main(args,"",output=False)
        impulselist[ind] = params['impulse']
        milist[ind] = params['mi']
        heightflist[ind] = params['heightfuel']
        heightolist[ind] = params['heightox']
        newheightlist[ind] = params['newheight']
        thrustlist[ind] = params['thrust']
        mdotlist[ind] = params['mdot']
        mdotoxlist[ind] = params['mdot_ox']
        mdotfuellist[ind] = params['mdot_fuel']
        isplist[ind] = params['isp']
        maxtwglist[ind] = params['twg_max']
        maxtwclist[ind] = params['twc_max']
        wstructlist[ind] = params['wstruct']
        lambdalist[ind] = params['L']
        totalmasslist[ind]=params['M']
        ind = ind+1


plt.plot(percents,impulselist,'g')
#plt.title("percents (psi) vs impulse (ns)")
plt.title("Ethanol Percent vs impulse (ns)")
plt.savefig(os.path.join(path, "impulse.png"))
plt.figure()
plt.plot(percents,isplist,'m')
plt.title("Ethanol Percent vs isp (s)")
plt.savefig(os.path.join(path, "isp.png"))
plt.figure()
plt.plot(percents,milist,'r')
plt.title("Ethanol Percent vs inert mass (kg)")
plt.savefig(os.path.join(path, "inertmass.png"))
plt.figure()
plt.plot(percents,thrustlist/const.lbToN)
plt.title("Ethanol Percent vs sea levelthrust(lb)")
plt.savefig(os.path.join(path, "thrust.png"))
plt.figure()
plt.plot(percents,newheightlist,'k')
plt.title("Ethanol Percent vs heights (m)")
plt.savefig(os.path.join(path, "height.png"))
plt.figure()
plt.plot(percents,mdotlist,'k', label='total')
plt.plot(percents,mdotoxlist,'b', label='ox')
plt.plot(percents,mdotfuellist,'r', label='fuel')
plt.title("Ethanol Percent vs mass flow rate (kg/s)")
plt.legend()
plt.savefig(os.path.join(path, "mdot.png"))
plt.figure()

plt.plot(percents,maxtwglist,'r', label='max twg')
plt.plot(percents,maxtwclist,'b', label='max twc')
plt.title("Ethanol Percent vs max wall temp (kelvin)")
plt.savefig(os.path.join(path, "chambertemps.png"))
plt.legend()

plt.show()
"""
"""
args = {
        'thrust': 4000 * const.lbToN,  # Newtons, INITIAL GUESS
        'time': 40,  # s, INIT GUESS
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        # 'phi':1,
        'TWR' : 4,
        'cr': 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa}
configtitle="NitroEthanol75_300psi_pressuredrop_TWR4"
path=os.path.join("Configs","Tradestudies",configtitle)
os.makedirs(path,exist_ok=False)

#twrs = np.arange(150,751,25)*const.psiToPa
dps = np.arange(50,300,25)
impulselist = np.zeros(dps.size)
milist = np.zeros(dps.size)
heightflist = np.zeros(dps.size)
heightolist = np.zeros(dps.size)
newheightlist=np.zeros(dps.size)
thrustlist = np.zeros(dps.size)
mdotlist = np.zeros(dps.size)
mdotoxlist = np.zeros(dps.size)
mdotfuellist = np.zeros(dps.size)
isplist = np.zeros(dps.size)
maxtwglist = np.zeros(dps.size)
maxtwclist = np.zeros(dps.size) #HAD TO COMMENT OUT LINE 317 and 318 incoling system!!
wstructlist = np.zeros(dps.size)
lambdalist = np.zeros(dps.size)
totalmasslist = np.zeros(dps.size)
ind = 0
for perc in dps:
        args['dp']=perc* const.psiToPa 
        print(perc)
        params=main(args,"",output=False)
        impulselist[ind] = params['impulse']
        milist[ind] = params['mi']
        heightflist[ind] = params['heightfuel']
        heightolist[ind] = params['heightox']
        newheightlist[ind] = params['newheight']
        thrustlist[ind] = params['thrust']
        mdotlist[ind] = params['mdot']
        mdotoxlist[ind] = params['mdot_ox']
        mdotfuellist[ind] = params['mdot_fuel']
        isplist[ind] = params['isp']
        maxtwglist[ind] = params['twg_max']
        maxtwclist[ind] = params['twc_max']
        wstructlist[ind] = params['wstruct']
        lambdalist[ind] = params['L']
        totalmasslist[ind]=params['M']
        ind = ind+1


plt.plot(dps,impulselist,'g')
#plt.title("dps (psi) vs impulse (ns)")
plt.title("Pressure  (psI) vs impulse (ns)")
plt.savefig(os.path.join(path, "impulse.png"))
plt.figure()
plt.plot(dps,isplist,'m')
plt.title("Pressure  (psI) vs isp (s)")
plt.savefig(os.path.join(path, "isp.png"))
plt.figure()
plt.plot(dps,milist,'r')
plt.title("Pressure  (psI) vs inert mass (kg)")
plt.savefig(os.path.join(path, "inertmass.png"))
plt.figure()
plt.plot(dps,thrustlist/const.lbToN)
plt.title("Pressure  (psI) vs sea levelthrust(lb)")
plt.savefig(os.path.join(path, "thrust.png"))
plt.figure()
plt.plot(dps,newheightlist,'k')
plt.title("Pressure  (psI) vs heights (m)")
plt.savefig(os.path.join(path, "height.png"))
plt.figure()
plt.plot(dps,mdotlist,'k', label='total')
plt.plot(dps,mdotoxlist,'b', label='ox')
plt.plot(dps,mdotfuellist,'r', label='fuel')
plt.title("Pressure  (psI) vs mass flow rate (kg/s)")
plt.legend()
plt.savefig(os.path.join(path, "mdot.png"))
plt.figure()

plt.plot(dps,maxtwglist,'r', label='max twg')
plt.plot(dps,maxtwclist,'b', label='max twc')
plt.title("Pressure  (psI) vs max wall temp (kelvin)")
plt.savefig(os.path.join(path, "chambertemps.png"))
plt.legend()

plt.show()
"""
args = {
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 10 * const.psiToPa,
       # 'phi':1,
        'cr': None,
        'TWR' : 4,
        'lstar': 1,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 / 2,
        'dp': 150 * const.psiToPa,
        'impulseguess' :  493555.8,
        'rc' : .11} #623919}
configtitle="TEST"
path=os.path.join("Configs","Tradestudies",configtitle)
os.makedirs(path,exist_ok=False)

#twrs = np.arange(150,751,25)*const.psiToPa
pes = np.arange(8,16,.25)
impulselist = np.zeros(pes.size)
milist = np.zeros(pes.size)
heightflist = np.zeros(pes.size)
heightolist = np.zeros(pes.size)
newheightlist=np.zeros(pes.size)
thrustlist = np.zeros(pes.size)
mdotlist = np.zeros(pes.size)
mdotoxlist = np.zeros(pes.size)
mdotfuellist = np.zeros(pes.size)
isplist = np.zeros(pes.size)
maxtwglist = np.zeros(pes.size)
maxtwclist = np.zeros(pes.size) #HAD TO COMMENT OUT LINE 317 and 318 incoling system!!
wstructlist = np.zeros(pes.size)
lambdalist = np.zeros(pes.size)
totalmasslist = np.zeros(pes.size)
ind = 0
for perc in pes:
        args['pe']=perc* const.psiToPa 
        print(perc)
        params=main(args,"",output=False)
        impulselist[ind] = params['impulse']
        milist[ind] = params['mi']
        heightflist[ind] = params['heightfuel']
        heightolist[ind] = params['heightox']
        newheightlist[ind] = params['newheight']
        thrustlist[ind] = params['thrust']
        mdotlist[ind] = params['mdot']
        mdotoxlist[ind] = params['mdot_ox']
        mdotfuellist[ind] = params['mdot_fuel']
        isplist[ind] = params['isp']
        maxtwglist[ind] = params['twg_max']
        maxtwclist[ind] = params['twc_max']
        wstructlist[ind] = params['wstruct']
        lambdalist[ind] = params['L']
        totalmasslist[ind]=params['M']
        ind = ind+1


plt.plot(pes,impulselist,'g')
#plt.title("pes (psi) vs impulse (ns)")
plt.title("Pressure  (psI) vs impulse (ns)")
plt.savefig(os.path.join(path, "impulse.png"))
plt.figure()
plt.plot(pes,isplist,'m')
plt.title("Pressure  (psI) vs isp (s)")
plt.savefig(os.path.join(path, "isp.png"))
plt.figure()
plt.plot(pes,milist,'r')
plt.title("Pressure  (psI) vs inert mass (kg)")
plt.savefig(os.path.join(path, "inertmass.png"))
plt.figure()
plt.plot(pes,thrustlist/const.lbToN)
plt.title("Pressure  (psI) vs sea levelthrust(lb)")
plt.savefig(os.path.join(path, "thrust.png"))
plt.figure()
plt.plot(pes,newheightlist,'k')
plt.title("Pressure  (psI) vs heights (m)")
plt.savefig(os.path.join(path, "height.png"))
plt.figure()
plt.plot(pes,mdotlist,'k', label='total')
plt.plot(pes,mdotoxlist,'b', label='ox')
plt.plot(pes,mdotfuellist,'r', label='fuel')
plt.title("Pressure  (psI) vs mass flow rate (kg/s)")
plt.legend()
plt.savefig(os.path.join(path, "mdot.png"))
plt.figure()

plt.plot(pes,maxtwglist,'r', label='max twg')
plt.plot(pes,maxtwclist,'b', label='max twc')
plt.title("Pressure  (psI) vs max wall temp (kelvin)")
plt.savefig(os.path.join(path, "chambertemps.png"))
plt.legend()

plt.show()