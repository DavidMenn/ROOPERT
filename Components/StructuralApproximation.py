"""Here will hopefully lie the code that takes overall
rocket params and turns them into weights! Currently just copypasted what I have from matlab
and full of errors and over-simplifications
when you rewrite, fill out mass_approx_new and then rename it to mass_approx and
change mass_approx to mass_approx_OLD"""
import sys
sys.path.insert(1,"./")
import math
from Toolbox.Constant import psiToPa, lbToKg
import os
def mass_approx(Pc,dp, OD, rho_f,rho_ox, thrust, isp, burntime, rm, printoutput=False, outputdir = ""):
    strcutresweight = 50.6*lbToKg;
    propweight =89.72*lbToKg;
    cfuelweight =5.3*lbToKg;
    coxweight=16.7 *lbToKg;
    cpresweight=14.5*lbToKg;
    cpressys= 12.21*lbToKg;
    cfoxhole=7.7*lbToKg;
    cav=15.03*lbToKg;
    cTC=16.66 *lbToKg;
    currentheight = 253.47*.0254;
    fuelheight=8.72*.0254;
    oxheightcurrent=37.824*.0254;
    presheightcurrent=16*.0254;

    diam_tank	= (OD-.1
                    )*.0254;
    ulox	 =0.05;
    ulfuel	=0.05;
    ulpres	=0.1;
    mf = thrust/isp/9.81; #just pass mass flow rate
    mf_ox=(rm*mf/(rm+1));
    mf_fuel = mf/(rm+1);
    vol_ox	= mf_ox*burntime/rho_ox;
    vol_fuel	= mf_fuel*burntime/rho_f;
    P_pres =	4500*psiToPa;
    P_tank =	Pc+dp;
    vol_pres =P_tank*(vol_fuel+vol_ox)/(P_pres-P_tank); #compresasbility


    heightox	=(vol_ox/(1-ulox))/(math.pi*(diam_tank/2)**2);
    heightfuel	=(vol_fuel/(1-ulfuel))/(math.pi*(diam_tank/2)**2);
    heightpres	=(vol_pres/(1-ulpres))/(math.pi*(diam_tank/2)**2);

    Sy = 40000*psiToPa;
    Fs	=2;
    rho_tank =	2710 #kg/m3
    t_prop	=P_tank*diam_tank/(2*Sy)*Fs;
    t_pres	=P_pres*diam_tank/(2*Sy)*Fs;

    weight_f	=((math.pi*diam_tank*heightfuel)+2*math.pi*(diam_tank/2)**2)*t_prop*rho_tank;
    weight_ox	=((math.pi*diam_tank*heightox)+2*math.pi*(diam_tank/2)**2)*t_prop*rho_tank;
    weight_pres	=((math.pi*diam_tank*heightpres)+2*math.pi*(diam_tank/2)**2)*t_pres*rho_tank;

    massFuelTank = 0
    massOxTank = 0
    while abs(weight_f-massFuelTank)>.01:
        weight_f=massFuelTank
        newheight=(currentheight-(fuelheight+oxheightcurrent+presheightcurrent))*8/OD+(heightox+heightfuel+heightpres);#(currentheight-(fuelheight+oxheightcurrent+presheightcurrent))*OD/8+(heightox+heightfuel+heightpres);
        Saratio=(newheight*math.pi/4*OD**2)/(currentheight*math.pi/4*8**2); #GET RID OF THIS
        mdotratio=mf/(4.72*0.453592);
        Pratio_prop=P_tank/(950*psiToPa);
        Pratio_pres=P_pres/(4500*psiToPa);

        global wtc, wav, wfox, wpresys, wpres, whelium, wstruct

        wtc=mdotratio*(cTC);
        wav=mdotratio*Pratio_prop*(cav); # FIX ALL THESE
        wfox=mdotratio*Pratio_prop*(cfoxhole);
        wpresys=mdotratio*Pratio_pres*cpressys;
        wpres=(vol_pres/0.0055217851)*cpresweight; #make sure this is in metric
        whelium=vol_pres*44.472*P_pres/(31.02*10**6); # this is denisty from refprop for 4500 psi

        # # Ideally You would use these not an overall structures weight
        massNoseCone =  11.83*lbToKg # THIS IS A GUESS
        massFins = 4 # THIS IS ALSO A GUESS
        massRecovery = 10 # THIS IS ALSO A GUESS
        massRecoveryCoupler = 8*lbToKg;
        massPayloadCoupler = 2.49;
        massAirFrames = 19.962*newheight/currentheight #air frames mass scales with height!?
        massPayload = 2.5
        wstruct=massPayload+massAirFrames+massPayloadCoupler+massRecovery+massRecoveryCoupler+massFins+massNoseCone
        # massAVCoupler = 0;
        # massParachuteShroud = 0;
        # massAvionicsShroud = 0;
        # massAvionicsShroud = 0;
        # massIntertankBay1Shroud = 0;

        #wstruct=Saratio*strcutresweight#1.369*Saratio*strcutresweight; #The factor for structural weight is in the error in newheight estimate lol its like 1.2


        # BREAK DOWN STUCTURES WEIGHT INTO COMPONENTS
        oxTankThickness,fuelTankThickness, heightfuel, heightox, massOxTank, \
        massFuelTank,fuelTankLengthTotal,oxTankLengthTotal = MetalTankMasses((OD)*0.0254, P_tank, weight_ox, vol_ox,vol_fuel, mf_ox*burntime)

    weight_f=massFuelTank
    weight_ox = massOxTank



    mis=wtc+wav+wfox+wpresys+wpres+wstruct+weight_f+weight_ox+whelium;
    fuelmass=mf_fuel*burntime;
    oxmass=mf_ox*burntime;
    lambdas=(fuelmass+oxmass)/(mis+fuelmass+oxmass);
    totalmasses=fuelmass+oxmass+mis;
    if printoutput:
        os.makedirs(outputdir, exist_ok=True)
        with open(os.path.join(outputdir,"mass_output.txt"),"a")  as f:
            print(f"vol_ox, {vol_ox}", file=f)
            print(f"vol_fuel, {vol_fuel}", file=f)
            print(f"P_pres, {P_pres}", file=f)
            print(f"P_tank, {P_tank}", file=f)
            print(f"vol_pres, {vol_pres}", file=f)
            print(f"t_prop: ox tank, {oxTankThickness}", file=f)
            print(f"weight_fueltank, {weight_f}", file=f)
            print(f"weight_oxtank, {weight_ox}", file=f)
            print(f"heightrocket, {newheight}", file=f)
            print(f"heightfuel, {heightfuel}", file=f)
            print(f"heightox, {heightox}", file=f)
            print(f"thrust chamber weight, {wtc}", file=f)
            print(f"avbay weight, {wav}", file=f)
            print(f"foxhole weight, {wfox}" , file=f)
            print(f"massNoseCone, {massNoseCone}", file=f)
            print(f"massFins, {massFins}", file=f)
            print(f"massRecovery, {massRecovery}", file=f)
            print(f"massRecoveryCoupler, {massRecoveryCoupler}", file=f)
            print(f"massPayloadCoupler, {massPayloadCoupler}", file=f)
            print(f"massAirFrames, {massAirFrames}", file=f)
            print(f"massPayload, {massPayload}", file=f)
            print( f"Structures weight, {wstruct}", file=f)
            print(f"inert mass, {mis}", file=f)
            print(f"fuelmass, {fuelmass}", file=f)
            print(f"oxmass, {oxmass}", file=f)
            print(f"lambda (massfrac), {lambdas}", file=f)
            print(f"totalmass, {totalmasses}", file=f)


    return mis, lambdas, totalmasses, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank

def mass_approx_NEW():
    """Feel free to use whatever you want as inputs, though try to keep it limited
    and to things that can be calculated from first order calcs
    inputs : chamber pressure, thrust, burntime, propellant props, whatever else
    outputs : dry mass, mass fraction, total mass
    Note that in the matlab version i had this working with matrices to do optimization faster
    Idk if thats something we want to pursue in python so you can hjust reutrn doubles"""



#function [oxTankThickness,fuelTankThickness, fuelLength, oxLength, massOxTank, massFuelTank,fuelTankLengthTotal,oxTankLengthTotal] = MetalTankMasses(OD,tP,massOxTankGuess,oxVol,fuelVol)
def MetalTankMasses(OD,tP_PA,massOxTankGuess,oxVol,fuelVol, oxWeight):
    """
    Inputs:
    OD = Tank OD [m]
    tP_PA = Tank Pressure [pa]
    oxVol = volume for ox tank
    fuelVol = volume for fuel tank
    """

    # Calculating Load Conditions
    # Pressure
    tankOuterRadius = OD/2; # Tank OR from OD
    YieldStrength_Aluminum6061 = 241*10**6; # PA (different for cryo, just using normal for now)
    FOS = 2;    # Based on IREC Requirnments for COPV's - Subject to change
    maxStress = YieldStrength_Aluminum6061/FOS;
    maxGs = 7;

    # Axial Loading
    # Mass Assumptions (kg)
    massAvBayPlumbing = wav
    massFoxHolePlumbing = wfox
    massPressTank = wpres
    massPressSystem = wpresys
    massHelium = whelium
    massStructures = wstruct # whole rocket - ideally you would only consider structures above each tank

    #fuelWeight = 102;
    #oxWeight = 213.45;

    # # Ideally You would use these not an overall structures weight
    # massNoseCone = 0;
    # massFins = 0;
    # massRecovery = 0;
    # massRecoveryCoupler = 0;
    # massPayloadCoupler = 2.49;
    # massAVCoupler = 0;
    # massParachuteShroud = 0;
    # massAvionicsShroud = 0;
    # massAvionicsShroud = 0;
    # massIntertankBay1Shroud = 0;
    # massOverFuel = massNoseCone + massPayloadCoupler + massParachuteShroud + ...
    # massRecoveryCoupler + massAvBayPlumbing + massAvionicsShroud + massPressTank+ massHelium + massPressSystem;
    # massOverOx = massOverFuel + massFoxHolePlumbing + massIntertankBay1Shroud;

    massOverFuel = massStructures + massAvBayPlumbing + massPressTank+ massHelium + massPressSystem ;
    massOverOx = massOverFuel + massFoxHolePlumbing + oxWeight + massOxTankGuess;

    # Axial Load From mass estimates
    OxAxialLoad = massOverOx*maxGs*9.8; # N
    FuelAxialLoad = massOverFuel*maxGs*9.8; # N

    # Calculate Tank Thickness
    oxTankThickness = thinWallTankThicknessApproximation(maxStress,OxAxialLoad,tP_PA,tankOuterRadius);
    fuelTankThickness = thinWallTankThicknessApproximation(maxStress,FuelAxialLoad,tP_PA,tankOuterRadius);

    # Find Tank IR
    fuelTankIR = tankOuterRadius-fuelTankThickness;
    oxTankIR = tankOuterRadius-oxTankThickness;

    # Internal Length From Volume
    oxLength = oxVol/(math.pi*(oxTankIR)**2);
    fuelLength = fuelVol/(math.pi*(fuelTankIR)**2);

    # Tank Mass
    massBulkhead = 6.57; # kg, i have to believe this
    massOxTube = massCylinder(oxTankIR,tankOuterRadius,oxLength);
    massFuelTube = massCylinder(fuelTankIR,tankOuterRadius,fuelLength);
    massOxTank = 2*massBulkhead+massOxTube;
    massFuelTank = 2*massBulkhead+massFuelTube;

    # Tank Heights with bulkhead (approximate)
    bulkheadLengthAboveTube = .17;  # m
    oxTankLengthTotal = bulkheadLengthAboveTube*2 + oxLength;
    fuelTankLengthTotal = bulkheadLengthAboveTube*2 + fuelLength;
    return oxTankThickness,fuelTankThickness, fuelLength, oxLength, massOxTank, massFuelTank,fuelTankLengthTotal,oxTankLengthTotal
    # Uniform internal pressure q (ends capped) + Axial Load
    # Inputs: maximum yield strength of material (includes factor of safety),
    # Max Operating pressure, estimated axial load, tank radius
    # Outputs: Approximate tank wall thickness
def thinWallTankThicknessApproximation(maxSigma,AxialLoad,MEOP,R):
    q = MEOP;
    p = AxialLoad;
    tAxial = (1/maxSigma) * (q*R/2+p);
    tHoop = q*R/maxSigma;
    thickness = max([tAxial,tHoop]); #Does this make sense? why arent we doing vector sum
    return thickness
# Inputs: Inner Radius and Outer Radius of tank
def massCylinder(IR,OR,length):
    density = 2720; # kg/m**3
    Area = math.pi*(OR**2-IR**2);
    mass = Area*length*density; # kg
    return mass
