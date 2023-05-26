import unittest
from Analysis import FirstOrderCalcs as FAC
from Toolbox import Constant as const
import math

class TestFOC(unittest.TestCase):
    def setUp(self) -> None:
        self.args = {
        'thrust': 1200 * const.lbToN,  # Newtons
        'time': 7.5,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 10 * const.psiToPa,
        'phi':1,
       # 'cr': 5.295390217,
        'lstar': 1.24,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 *2,
        'dp': 150 * const.psiToPa,
        'impulseguess' :  495555.24828424345,
        'rc' : .08,
        'thetac' : (35*math.pi/180),
        'isp_efficiency' : .9}
        self.params  = FAC.SpreadsheetSolver(self.args)
    def test_cr(self):
        tempargs = self.args.copy()
        tempargs['cr'] = self.params['cr']
        tempargs['rc'] = None
        tempparams = FAC.SpreadsheetSolver(tempargs)
        self.assertAlmostEqual(tempparams['rc'],self.args['rc'])

