"""This is the essence of Phase 1.
The idea is you pass it a dicitonary of parameters you know, and it'll fill out the dictionary
Look at the dictionary paramnames to see what variable options are supported
The variable calculations are stored in the functions after SpreadsheetSovler
An example of how to use this module:

import FirstOrderCalcs as FAC
args = {
        'thrust': 5000 * const.lbToN,  # Newtons
        'time': 30,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 350 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
        'phi': 1,
        'fuelname': 'Ethanol',
        'oxname': 'LOX',
        'throat_radius_curvature': .02}

    params = FAC.SpreadsheetSolver(args)
"""

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
def SpreadsheetSolver(args : dict, defaults : dict = None):
    paramnames = [
    'thrust', #Newtons
    'time', #s
    'impulse', #N*s
    'rho_ox', #Kg/M^3
    'rho_fuel', #Kg/M^3
    'pc', #Pa, if cr is specified, this is Pressure at end of combustor
    'pinj',#Pa, only useful if you specify CR, otherwise assumed to be pc
    'pe', #Pa
    'g', #m/s^2
    'rm', #o/f by weight
    'phi', #ratio from stoich (1 is stoich, >1 is fuel rich)
    'at', # m^2, area of throat
    'rt', # m, radius of throat
    'cr', # contraction ratio
    'rc', # m, combustion chamber radius
    'ac', # m^2, area combustion chamber
    'l_star', # m, volume cc/area throat
    'mol_weight', # kg/mol
    'gamma', # in cc
    'gamma_exit', # in exit
    'gamma_throat', # in throat
    'isp', # s
    'temp_c', # K, chamber temp
    'rg', # specific gas constant (SI units if what they are)
    'pr_throat',
    'rho_throat',
    'temp_e',
    'v_exit',
    'a_exit',
    'mach_exit',
    'temp_throat',
    'p_throat',
    'v_throat',
    'mdot',
    'mdot_ox',
    'mdot_fuel',
    'er',
    'cstar',
    'cf',
    'c_eff',
    'rho_av',
    'vc',
    'theta_con',
    'lc',
    'theta_div',
    'ln_conical',
    'ln_bell',
    'throat_radius_curvature',
    'ae',
    're',
    'nv',
    'nvstar',
    'nf',
    'nw',
    'fuelname',
    'oxname',
    'CEA',
    'pambient',
    'cf_efficiency', # Huzel and Huang page 16
    'isp_efficiency']

    params = dict.fromkeys(paramnames)

    params['g'] = 9.81

    #HERE IS WHERE YOU SET DEFAULT VALUES EARLY SO THINGS DONT BREAK IF YOU FORGET TO PASS STUFF
    if params['pambient'] is None:
        params['pambient']=14.7*const.psiToPa
    if params['pe'] is None:
        params['pe'] = params['pambient'] * 6894.76
    if params['cf_efficiency'] is None:
        params['cf_efficiency'] = .95
    if params['isp_efficiency'] is None:
        params['isp_efficiency']= .95

    for arg in list(args):
        try:
            params[arg] = args[arg]
        except:
            try:
                print("Parameter " + arg + " isn't supported or is mispelled. Did you mean " +
                      difflib.get_close_matches(arg, paramnames)[0] + "? That's what I'm going to use!")
                params[difflib.get_close_matches(arg, paramnames)[0]] = args[arg]
            except:
                print("Parameter" + arg + "isn't supported or is mispelled")

    if params['CEA'] is None: # need to have cea object to do other stuff sio make sure to get this done first
        try:
            RA.makeEthanolBlend(int(regex.search(r'\d+', params['fuelname']).group()))
        except:
            print(params['fuelname'])
        if params['cr'] is None:
            params['CEA'] = CEA_Obj(oxName=params['oxname'], fuelName=params['fuelname'], useFastLookup=1, isp_units='sec', cstar_units='m/s',
                                    pressure_units='Pa', temperature_units='K', sonic_velocity_units='m/s',
                                    enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units='J/kg-K',
                                    viscosity_units='millipoise', thermal_cond_units='W/cm-degC', fac_CR=None,
                                    make_debug_prints=False)
        else:
            params['CEA'] = CEA_Obj(oxName=params['oxname'], fuelName=params['fuelname'],useFastLookup=1,  isp_units='sec',
                                    cstar_units='m/s',
                                    pressure_units='Pa', temperature_units='K', sonic_velocity_units='m/s',
                                    enthalpy_units='J/kg', density_units='kg/m^3', specific_heat_units='J/kg-K',
                                    viscosity_units='millipoise', thermal_cond_units='W/cm-degC', fac_CR=params['cr'],
                                    make_debug_prints=False)
            PinjOverPcombEstimate = 1.0 + 0.54 / params['cr'] ** 2.2 # FROM https://rocketcea.readthedocs.io/en/latest/finite_area_comb.html
            params['pinj']=params['pc']*params['CEA'].cea_obj.get_Pinj_over_Pcomb(
                Pc=params['CEA'].Pc_U.uval_to_dval(params['pc'])*PinjOverPcombEstimate, MR=1.0, fac_CR=params['cr'])

        if params['rm'] is None:
            params['rm'] = rm(params)
        if params['pc'] is None:
            params['pc'] = pc(params)
        if params['er'] is None:
            params['er'] = er(params)


        molandgam = params['CEA'].get_Chamber_MolWt_gamma(Pc=params['pc'], MR=params['rm'], eps=params['er']),  # kg/mol
        params['mol_weight'] = molandgam[0][0]
        params['gamma'] = molandgam[0][1]
        molandgam = params['CEA'].get_exit_MolWt_gamma(Pc=params['pc'], MR=params['rm'], eps=params['er']),  # kg/mol
        params['gamma_exit'] = molandgam[0][1]  # in exit
        molandgam = params['CEA'].get_Throat_MolWt_gamma(Pc=params['pc'], MR=params['rm'], eps=params['er']),  # kg/mol
        params['gamma_throat'] = molandgam[0][1]
        temps = params['CEA'].get_Temperatures(Pc=params['pc'], MR=params['rm'], eps=params['er'], frozen=0,
                                               frozenAtThroat=0)
        params['isp'] = params['isp_efficiency']*\
        params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'], Pamb=params['pambient'],
                                           frozen=0, frozenAtThroat=0)[0]
        params['temp_c'] = temps[0]
        throatTrasnport = params['CEA'].get_Throat_Transport(Pc=params['pc'], MR=params['rm'], eps=params['er'],
                                                             frozen=0)
        params['rg'] = const.R/ params['mol_weight'],  # specific gas constant (SI units if what they are)
        params['pr_throat'] = throatTrasnport[3]
        params['rho_throat'] = \
        params['CEA'].get_Densities(Pc=params['pc'], MR=params['rm'], eps=params['er'], frozen=0, frozenAtThroat=0)[1]
        params['temp_e'] = temps[2]
        sonicvelos = params['CEA'].get_SonicVelocities(Pc=params['pc'], MR=params['rm'], eps=params['er'], frozen=0,
                                                       frozenAtThroat=0)
        params['mach_exit'] = params['CEA'].get_MachNumber(Pc=params['pc'], MR=params['rm'], eps=params['er'], frozen=0,
                                                           frozenAtThroat=0)

        params['a_exit'] = sonicvelos[2]
        params['v_exit'] = params['mach_exit'] * params['a_exit']

        params['temp_throat'] = temps[1]
        params['p_throat'] = params['pc']*params['CEA'].get_Throat_PcOvPe(Pc=params['pc'], MR=params['rm'])
        params['v_throat'] = sonicvelos[1]

        chambertransport = params['CEA'].get_Chamber_Transport(Pc=params['pc'], MR=params['rm'], eps=params['er'], frozen=0)
        params['prns'] = chambertransport[3]
        params['viscosityns'] = chambertransport[1]*10**(-6)
        params['cpns'] = chambertransport[0]

        params['cstar'] = params['isp_efficiency']*params['CEA'].get_Cstar(Pc=params['pc'], MR=params['rm'])

    for param in list(params):
        if params[param] is None:
            try:
                params[param] = globals()[param](params)
            except:
                print('Couldnt calculate ' + param + " first time around")

    #Do this twice so that if we need calcualted values for other values we can get them the second time around
    for param in list(params):
        if params[param] is None:
            try:
                params[param] = globals()[param](params)
            except:
                pass
    for param in list(params):
        if params[param] is None:
            try:
                params[param] = globals()[param](params)
            except:
                pass
    for param in list(params):
        if params[param] is None:
            try:
                params[param] = globals()[param](params)
            except:
                pass
    for param in list(params):
        if params[param] is None:
            try:
                params[param] = globals()[param](params)
            except:
                pass
    for param in list(params):
        if params[param] is None:
            try:
                params[param] = globals()[param](params)
            except:
                print(param + " not yet supported")
    return params

def thrust(params):
    try:
        return params['impulse']/params['time']
    except:
        print('Could not calculate Thrust')
        return None
def time(params): #s
    try:
        return params['impulse']/params['thrust']
    except:
        print('Could not calculate time')
        return None
def impulse(params): #N*s
    try:
        return params['thrust']*params['time']
    except:
        print('Could not calculate impulse')
        return None
def rho_ox(params): #Kg/M^3 IMPLEMENT THIS USING FLUID PROPS
    try:
        pObj = get_prop(params['oxname'])
        return 1000*pObj.SG_compressed(pObj.TdegRAtPsat(40),
                                  params['pc'] / const.psiToPa)
    except:
        print("rocketProps busted lol")
    return None
def rho_fuel(params): #Kg/M^3MPLEMENT THIS USING FLUID PROPS
    try:
        pObj = get_prop(params['fuelname'])
        return 1000 * pObj.SG_compressed(293*const.degKtoR,
                                         params['pc'] / const.psiToPa)
    except AttributeError:
        print('Fuel does not exist, assuming its a water ethanol blend')
        ethpercent = int(regex.search(r'\d+', params['fuelname']).group())/100
        pObj = get_prop("ethanol")
        ethdens = 1000 * pObj.SG_compressed(293*const.degKtoR,
                                         params['pc'] / const.psiToPa)
        pObj = get_prop("water")
        waterdens = 1000 * pObj.SG_compressed(293 * const.degKtoR,
                                            params['pc'] / const.psiToPa)
        return ethpercent*ethdens +(1-ethpercent)*waterdens
    except:
        print("rocketProps busted lol")
    return None
def pc(params): #Pa, calculate using throat conditions and backtrack
    print('You should input PC, will set up calculator later, now its just 700 psi (4.82 mpa)')
    return 700*6894.76
def pe(params): #Pa, this hsould never be called
    print('this should not be called this is a default value')
    return 14.7
def g(params): #m/s^2 his hsould never be called
    print('this should not be called this is a default value')
def rm(params): #o/f by weight
    try:
        return params['CEA'].getMRforER(ERphi=params['phi'])
    except:
        print('INSTALL ROCKETCEA IDIOT')



def phi(params): #ratio from stoich (1 is stoich, >1 is fuel rich)
    temp = params['CEA'].get_eqratio(Pc=params['pc'], MR=params['rm'])
    return temp[0]
def at(params): # m^2, area of throat assume choked at throat
    #return IE.AreaForChokedFlow(params['p_throat'],params['temp_throat'],params['gamma_throat'],params['mdot'], const.R/params['mol_weight'])
    if params['cr'] is None:
        #return IE.AreaForChokedFlow( IE.totalP(params['p_throat'],params['gamma_throat'],1), IE.totalT(params['temp_throat'],params['gamma_throat'],1) ,params['gamma_throat'], params['mdot'],
        #                        const.R / params['mol_weight']*1000)
        return params['thrust'] / params['cf'] / params['pc']
    else:
        #return IE.AreaForChokedFlow( IE.totalP(params['pc'],params['gamma'],IE.machFromP(params['pinj'],params['pc'],params['gamma'])),
        #                            IE.totalT(params['temp_c'],params['gamma'],IE.machFromP(params['pinj'],params['pc'],params['gamma'])) ,params['gamma_throat'], params['mdot'],
        #                      const.R / params['mol_weight']*1000)
        return params['thrust']/params['cf']/params['pc']
def rt(params): # m, radius of throat
    return math.sqrt(params['at']/const.pi)
"""THESE NEXT THREE ARE SUPER PRELIMINARY SIZINGS BASED OFF OF DIAGRAMS IN HUZEL AND HUANG PAGE 73"""
def cr(params): # contraction ratio
    CRInterpolator = interpolate.interp1d([i*const.inToM for i in [.2, .5, 1, 2, 5]],
                                          [i*const.inToM for i in [15, 9, 5.5, 4, 3]], kind='linear') #CR vs Throat Diam, visually copied from textbook
    print(CRInterpolator(params['rt']*2))
    return CRInterpolator(params['rt']*2)
def lc(params):
    CLInterpolator = interpolate.interp1d([i * const.inToM for i in [.2, .5, 1, 2, 5, 10]],
                                          [i * const.inToM for i in [1.8,3,4,5,7.5, 10]],
                                          kind='linear')  # CL vs Throat Diam, visually copied from textbook
    return CLInterpolator(params['rt'] * 2)
def rc(params): # m, combustion chamber radius
    if params['cr'] is None:
        return math.sqrt(params['l_star']*params['at']/(params['lc']*const.pi))
    else:
        return math.sqrt(params['at']*params['cr']/const.pi)
def ac(params): # m^2, area combustion chamber
    if params['cr'] is None:
        return const.pi*(params['l_star'] * params['at'] / (params['lc']))
    else:
        return params['at']*params['cr']
def l_star(params): # m, volume cc/area throat
    return 40*const.inToM # THIS IS KERALOX LSTAR FROM HUZEL AND HUANG PAGE 72 WE HSOULD DO MORE RESEARCH
def mdot(params):
    return params['thrust']/params['isp']/(const.gravity)
def mdot_ox(params):
    return params['rm']*params['mdot']/(params['rm']+1)
def mdot_fuel(params):
    return params['mdot']/(params['rm']+1)
def er(params):
    return params['CEA'].get_eps_at_PcOvPe(Pc=params['pc'], MR=params['rm'], PcOvPe=params['pc']/params['pambient'], frozen=0, frozenAtThroat=0)

def cf(params):
    return params['cf_efficiency']*((2*params['gamma']**2)/(params['gamma']-1)*(2/(params['gamma']+1))**((params['gamma']+1)/(params['gamma']-1))*(1-(params['pe']/params['pc'])**((params['gamma']-1)/params['gamma'])))**.5+params['er']*(params['pe']-params['pambient'])/params['pc']
"""def c_eff(params):
def rho_av(params):
def vc(params):
def theta_con(params):
def theta_div(params):

def ln_bell(params):
"""
def ln_conical(params):
    return (params['re']-params['rt'])/math.sin(15*math.pi/180)
def ae(params):
    return params['at']*params['er']
def re(params):
    return math.sqrt(params['at']*params['er']/const.pi)
"""
def nv(params):
def nvstar(params):
def nf(params):
def nw(params):
def fuelname(params):
def oxname(params):
"""






