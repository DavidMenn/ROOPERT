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
import FirstOrderCalcs as FAC
from scipy import interpolate

args = {
    'thust': 3000 * const.lbToN,  # Newtons
    'time': 50,  # s
    # 'rho_ox' : 1141, #Kg/M^3
    # 'rho_fuel' : 842,
    'pc': 350 * const.psiToPa,
    'pe': 14.7 * const.psiToPa,
    'phi': 1,
    'mdot': 18,
    'fuelname': 'Ethanol_75',
    'oxname': 'LOX'}

params = FAC.SpreadsheetSolver(args)


'''TO DO:
Get Tank volumes
- structural approx
- max accelerations
- sizes of engine components
- weights, shock force?
- apogee height
'''