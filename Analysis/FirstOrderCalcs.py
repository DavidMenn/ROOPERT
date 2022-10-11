import numpy as np
import sympy as sym
import difflib
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import Toolbox.RocketCEAAssister as RA
import re as regex

def SpreadsheetSolver(args : dict, defaults : dict = None):
    paramnames = [
    'thrust', #Newtons
    'time', #s
    'impulse', #N*s
    'rho_ox', #Kg/M^3
    'rho_fuel', #Kg/M^3
    'pc', #Pa
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
    'ae',
    're',
    'nv',
    'nvstar',
    'nf',
    'nw',
    'fuelname',
    'oxname',
    'CEA']

    params = dict.fromkeys(paramnames)

    if defaults is None:
        params['g'] = 9.81
        params['pe'] = 14.7 * 6894.76

    for arg in list(args):
        try:
            params[arg] = args[arg]
        except:
            try:
                print("Parameter" + arg + "isn't supported or is mispelled. Did you mean " + difflib.get_close_matches(
                    arg, paramnames) + "? That's what I'm going to use!")
                params[difflib.get_close_matches(arg, paramnames)] = args[difflib.get_close_matches(arg, paramnames)]
            except:
                print("Parameter" + arg + "isn't supported or is mispelled")

    if params['CEA'] is None: # need to have cea object to do other stuff sio make sure to get this done first
        RA.makeEthanolBlend(int(regex.search(r'\d+', args['fuelname']).group()))
        if params['rm'] is None:
            params['rm'] = rm(params)
        if params['pc'] is None:
            params['pc']  = pc(params)
        if params['er'] is None:
            params['er'] = er(params)
        params['CEA'] = CEA_Obj(oxName=args['oxname'], fuelName=args['fuelname'])
        (params['mol_weight'],params['gamma']) = params['CEA'].get_Chamber_MolWt_gamma(Pc=params['pc'], MR=params['rm'], eps=params['er']),  # kg/mol
        ,  # in cc
        params['gamma_exit'],  # in exit
        params['gamma_throat'],  # in throat
        params['isp'],  # s
        params['temp_c',  # K, chamber temp
        params['rg',  # specific gas constant (SI units if what they are)
        params['pr_throat',
        params['rho_throat',
        params['temp_e',
        params['v_exit',
        params['a_exit',
        params['mach_exit',
        params['temp_throat',
        params['p_throat',
        params['v_throat',


    for param in list(params):
        if params[param] is None:
            params[param] = locals()[param](params)

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
    print('Could not calculate RHO_OX')
    return None
def rho_fuel(params): #Kg/M^3MPLEMENT THIS USING FLUID PROPS
    print('Could not calculate RHO_FUel')
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
        return parms.CEA.getMRforER(ERphi=1)
    except:
        print('INSTALL ROCKETCEA IDIOT')



def phi(params): #ratio from stoich (1 is stoich, >1 is fuel rich)
    temp = get_eqratio(Pc=params.pc, MR=params.rm)
    return temp[0]
def at(params): # m^2, area of throat assume choked at throat
    return AreaForChokedFlow(params.CEA,T,gam,mdot,specificR)
def rt(params): # m, radius of throat
def cr(params): # contraction ratio
def rc(params): # m, combustion chamber radius
def ac(params): # m^2, area combustion chamber
def l_star(params): # m, volume cc/area throat
def mol_weight(params): # kg/mol
def gamma(params): # in cc
def gamme_exit(params): # in exit
def isp(params): # s
def temp_c(params): # K, chamber temp
def rg(params): # specific gas constant (SI units if what they are)
def pr_throat(params):
def rho_throat(params):
def temp_e(params):
def v_exit(params):
def a_exit(params):
def mach_exit(params):
def temp_throat(params):
def p_throat(params):
def v_throat(params):
def mdot(params):
def mdot_ox(params):
def mdot_fuel(params):
def er(params):
    er = 1
    temp = params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=er, Pamb=params['pe'], frozen=0, frozenAtThroat=0)
    while temp[1] == "UnderExpanded":
        er=er++
        temp = params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=er, Pamb=params['pe'], frozen=0,
                                                  frozenAtThroat=0)
    return er

def cstar(params):
def cf(params):
def c_eff(params):
def rho_av(params):
def vc(params):
def theta_con(params):
def lc(params):
def theta_div(params):
def ln_conical(params):
def ln_bell(params):
def ae(params):
def re(params):
def nv(params):
def nvstar(params):
def nf(params):
def nw(params):
def fuelname(params):
def oxname(params):







