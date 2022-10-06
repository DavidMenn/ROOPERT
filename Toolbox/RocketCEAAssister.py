"""Purpose of this is to help format CEA output, make fuel blends, etc. prevent repetitive stuff
ASSUMES YOU ALREADY HAVE CREATED A CEAOBJ USING THE FOLLOWING: (everything from https://rocketcea.readthedocs.io/en/latest/)
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
C = CEA_Obj( oxName='LOX', fuelName='LH2')

import numpy as np
"""
from rocketcea.cea_obj import add_new_fuel
def makeEthanolBlend(ethanolPercent):
    card_str = f"""
    fuel C2H5OH(L)  C 1.0   H 6.0   O 1.0    wt%={ethanolPercent}
    h,cal=-66326.482      t(k)=298.15   rho=.789
    fuel   H20(L) H 2.0 O 1.0   wt%={100-ethanolPercent}
    h,cal=-66326.482    t(k)=298.15   rho=1
    """
    add_new_fuel(f"Ethanol_{ethanolPercent}", card_str)
