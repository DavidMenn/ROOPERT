"""Here will hopefully lie the code that takes overall rocket params and turns them into weights! Currently just copypasted what I have from matlab"""
import math
from Toolbox.Constant import psiToPa, lbToKg
def mass_approx(Pc,dp, OD, rho_f,rho_ox, thrust, isp, burntime, rm):
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

    diam_tank	= (OD-2)*.0254;
    ulox	 =0.05;
    ulfuel	=0.05;
    ulpres	=0.1;
    mf = thrust/isp/9.81;
    mf_ox=(rm*mf/(rm+1));
    mf_fuel = mf/(rm+1);
    vol_ox	= mf_ox*burntime/rho_ox;
    vol_fuel	= mf_fuel*burntime/rho_f;
    P_pres =	4500*psiToPa;
    P_tank =	Pc+dp;
    vol_pres =P_tank*(vol_fuel+vol_ox)/(P_pres-P_tank);

    heightox	=(vol_ox/(1-ulox))/(math.pi*(diam_tank/2)**2);
    heightfuel	=(vol_fuel/(1-ulfuel))/(math.pi*(diam_tank/2)**2);
    heightpres	=(vol_pres/(1-ulpres))/(math.pi*(diam_tank/2)**2);

    Sy = 40000*psiToPa;
    Fs	=2;
    rho_tank =	2710 #kg/m3
    t_prop	=P_tank*diam_tank/(Sy)*Fs;
    t_pres	=P_pres*diam_tank/(Sy)*Fs;

    weight_f	=((math.pi*diam_tank*heightfuel)+2*math.pi*(diam_tank/2)**2)*t_prop*rho_tank;
    weight_ox	=((math.pi*diam_tank*heightox)+2*math.pi*(diam_tank/2)**2)*t_prop*rho_tank;
    weight_pres	=((math.pi*diam_tank*heightpres)+2*math.pi*(diam_tank/2)**2)*t_pres*rho_tank;

    newheight=(currentheight-(fuelheight+oxheightcurrent+presheightcurrent))*OD/8+(heightox+heightfuel+heightpres);#(currentheight-(fuelheight+oxheightcurrent+presheightcurrent))*8/OD+(heightox+heightfuel+heightpres);
    Saratio=(newheight*math.pi/4*OD**2)/(currentheight*math.pi/4*8**2);
    mdotratio=mf/4.72;
    Pratio_prop=P_tank/(950*psiToPa);
    Pratio_pres=P_pres/(4500*psiToPa);

    wtc=mdotratio*Pratio_prop*(cTC);
    wav=mdotratio*Pratio_prop*(cav);
    wfox=mdotratio*Pratio_prop*(cfoxhole);
    wpresys=mdotratio*Pratio_pres*cpressys;
    wpres=(vol_pres/0.0055217851)*cpresweight;
    whelium=(vol_pres*P_pres/101000);
    wstruct=Saratio*strcutresweight; #The factor for structural weight is in the error in newheight estimate lol its like 1.2

    mis=wtc+wav+wfox+wpresys+wpres+wstruct+weight_f+weight_ox+whelium;
    fuelmass=mf_fuel*burntime;
    oxmass=mf_ox*burntime;
    lambdas=(fuelmass+oxmass)/(mis+fuelmass+oxmass);
    totalmasses=fuelmass+oxmass+mis;

    return mis, lambdas, totalmasses
