# Foreward
# If Pressure Drop calculated is >10% of absolute inlet pressure,
# compressibility must be taken into account. This is in the form of a Y^2
# value explained on Crane 410 6-4

# Fluid Variables
    # Input Variables

import math
import Toolbox.Constant as const


def pressureDrop(OxDensityMetric = 1090.79390625913,     # kg/m^3
    OxMassFlowMetric = 4.09234712442909, # kg/s
    OxTubeDiam = 3/4, #in  
    FuelDensityMetric = 843.4,       # kg/m^3
    FuelMassFlowMetric = 3.28127,# kg/s
    FuelTubeDiam = 3/4     ,
    params = None           ):# in

    if params is None:
        params = {'pc' : 300*const.psiToPa}
        
    pi = math.pi # lol matlab
    OxLength = 2;                           # length of plumbing line in feet
    OxViscosity = 1015.24 / 1000;           # micropoise to centipoise

               # kg/s
   
    FuelLength = 10;                        # ft
    FuelViscosity = 1.095;                  # centipoise, approximated as 100% ethanol

    bendRatio = 3;                          # radius of bend / diameter of tube; 
                                            # r/d decreases with increased tube size. 3 is an estimate

        # Constants
    g = 32.174; # ft/s^2
    roughness = 0.015 / 304.8; # absolute roughness of stainless steel tubing, input mm, output ft

        # Calculated Variables
    OxDensity = OxDensityMetric/16.018; # lb/ft^3
    OxTubeDiamMetric = OxTubeDiam/39.37; #in to meters
    OxTubeAreaMetric = (pi/4)*(OxTubeDiamMetric**2); #m^2
    OxVelocityFlowMetric = OxMassFlowMetric/(OxDensityMetric * OxTubeAreaMetric); #m/s
    OxVelocityFlow = OxVelocityFlowMetric*3.281;
    OxReynold = 124*OxTubeDiam*OxVelocityFlow*OxDensity/OxViscosity;

    FuelDensity = FuelDensityMetric/16.018; # lb/ft^3
    FuelTubeDiamMetric = FuelTubeDiam/39.37; #in to meters
    FuelTubeAreaMetric = (pi/4)*(FuelTubeDiamMetric**2); #m^2
    FuelVelocityFlowMetric = FuelMassFlowMetric/(FuelDensityMetric * FuelTubeAreaMetric); #m/s
    FuelVelocityFlow = FuelVelocityFlowMetric*3.281;
    FuelReynold = 124*FuelTubeDiam*FuelVelocityFlow*FuelDensity/FuelViscosity;

    # Friction Coefficient f & K(Crane 410 Page X)
    Oxf = 0.25/(math.log(roughness/(3.7*OxTubeDiam/12) + 5.74/(OxReynold**0.9)))**2;
    OxfTurb = 0.25/(math.log(roughness/(3.7*OxTubeDiam/12)))**2;

    KOxElbow = 30*OxfTurb;
    KOxBend = 12*OxfTurb; # Tabularly depends on bend radius (See page 
    KOxBallValve = 3*OxfTurb;
    KOxCheckValve = 420*OxfTurb;
    KOxEnter = 0.78;
    KFuelExit = 1.0;

    KOx = KOxEnter + KOxBallValve + KOxCheckValve + KFuelExit;

    Fuelf = 0.25/(math.log(roughness/(3.7*OxTubeDiam/12) + 5.74/(FuelReynold**0.9)))**2;
    FuelfTurb = 0.25/(math.log(roughness/(3.7*FuelTubeDiam/12)))**2;

    KFuelElbow = 30*FuelfTurb;
    KFuelBend = 12*FuelfTurb;
    KFuelBallValve = 3*FuelfTurb;
    KFuelCheckValve = 420*FuelfTurb;
    KFuelEnter = 0.78;
    KFuelExit = 1.0;

    KFuel = KFuelEnter + 2*KFuelBend + KFuelBallValve + KFuelCheckValve + KFuelExit;

    # Pressure Drop (Crane 410 Page 6-4, Equation 6-22)
    OxPressureDropFriction = 0.001295*Oxf*OxLength*OxDensity*(OxVelocityFlow**2)/OxTubeDiam;
    OxPressureDropFittings = (OxDensity/144)*(KOx)*(OxVelocityFlow**2)/(2*g);

    OxPressureDropTotal = OxPressureDropFriction + OxPressureDropFittings

    FuelPressureDropFriction = 0.001295*Fuelf*FuelLength*FuelDensity*(FuelVelocityFlow**2)/FuelTubeDiam;
    FuelPressureDropFittings = (FuelDensity/144)*(KFuel)*(FuelVelocityFlow**2)/(2*g);

    FuelPressureDropTotal = FuelPressureDropFriction + FuelPressureDropFittings

    dparray =  [(params['pc']*.2 + 50*const.psiToPa + FuelPressureDropTotal)*1.5, (params['pc']*.2 + OxPressureDropTotal)*1.5]
    dparray = [150*const.psiToPa, 150*const.psiToPa]
    return OxPressureDropTotal, FuelPressureDropTotal, dparray
        # Note: pressureDropFriction is based on wall friction due to length of
        # straight pipe. pressureDropFittings is due to the turbulence caused
        # by the flow obstruction of internal geometries. It is not clear
        # whether the length of the fittings should be subtracted from
        # pressureDropFriction, or if the double-calculation is necessary.
