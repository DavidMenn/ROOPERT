"""Here will lie the wall temps. Early on we're just going to do throat temp just to get sizing down"""


# ns stands for nozzle start, use finite area combustor CEA to find it (or just use pc and tcomb)
#This equation neglects g because you dont need it for si units (i think)
def Bartz(D, viscosityns, prandtlns, cpns, pcns, cstar, throatRadiusCurvature, at, a, twg, tcns, mach, gamma):
    sigma = ((.5 * twg / tcns * (1 + (mach ** 2) * (gamma - 1) / 2) + .5) ** (-.68)) * (
                1 + (mach ** 2) * (gamma - 1) / 2) ** (-.12)
    return (.026 / (D ** .2) * (cpns * viscosityns ** .2) / (prandtlns ** .6) * (pcns / cstar) ** (.8) * (
                D / throatRadiusCurvature) ** .1) * (at / a) ** .9 * sigma
