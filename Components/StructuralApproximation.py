"""Here will hopefully lie the code that takes overall rocket params and turns them into weights! Currently just copypasted what I have from matlab"""

function [lambdas,mis,totalmasses,newheight] = mass_aprrox(Pc,dp, OD, rho_f,rho_ox, thrust, isp, burntime, rm)
    strcutresweight = 50.6; %structures
    propweight =89.72; % prop
    cfuelweight =5.3; %fuelt
    coxweight=16.7 ;%oxt
    cpresweight=14.5; %prestank
    cpressys= 12.21; %pressys
    cfoxhole=7.7; %foxholesys
    cav=15.03 ;%avsys
    cTC=16.66 ;%thrustchamber
    currentheight = 253.47; %curernt hweight
    fuelheight=8.72;
    oxheightcurrent=37.824;
    presheightcurrent=16;

    diam_tank	= OD-2;
    ulox	 =0.05;
    ulfuel	=0.05;
    ulpres	=0.1;
    mf = thrust./isp;
    mf_ox=(rm.*mf./(rm+1));
    mf_fuel = mf./(rm+1);
    vol_ox	= mf_ox.*burntime./rho_ox;
    vol_fuel	= mf_fuel.*burntime./rho_f;
    P_pres =	4500;
    P_tank =	Pc+dp;
    vol_pres =P_tank.*(vol_fuel+vol_ox)./(P_pres-P_tank);

    heightox	=(vol_ox./(1-ulox))./(pi.*(diam_tank./24)^2);
    heightfuel	=(vol_fuel./(1-ulfuel))./(pi.*(diam_tank./24)^2);
    heightpres	=(vol_pres./(1-ulpres))./(pi.*(diam_tank./24)^2);

    Sy = 40000;
    Fs	=2;
    rho_tank =	0.0975;
    t_prop	=P_tank.*diam_tank./(2.*Sy).*Fs;
    t_pres	=P_pres.*diam_tank./(2.*Sy).*Fs;

    weight_f	=((pi.*diam_tank.*heightfuel.*12)+2.*pi.*(diam_tank./2)^2).*t_prop.*rho_tank;
    weight_ox	=((pi.*diam_tank.*heightox.*12)+2.*pi.*(diam_tank./2)^2).*t_prop.*rho_tank;
    weight_pres	=((pi.*diam_tank.*heightpres.*12)+2.*pi.*(diam_tank./2)^2).*t_pres.*rho_tank;

    newheight=(currentheight-(fuelheight+oxheightcurrent+presheightcurrent)).*OD./8+(heightox+heightfuel+heightpres).*12;
    Saratio=(newheight.*pi./4.*OD^2)./(currentheight.*pi./4.*8^2);
    mdotratio=mf./4.72;
    Pratio_prop=P_tank./950;
    Pratio_pres=P_pres./4500;

    wtc=mdotratio.*Pratio_prop.*(cTC);
    wav=mdotratio.*Pratio_prop.*(cav);
    wfox=mdotratio.*Pratio_prop.*(cfoxhole);
    wpresys=mdotratio.*Pratio_pres.*cpressys;
    wpres=(vol_pres./0.195).*cpresweight;
    whelium=(vol_pres.*7.62);
    wstruct=Saratio.*strcutresweight*.8;

    mis=wtc+wav+wfox+wpresys+wpres+wstruct+weight_f+weight_ox+whelium;
    fuelmass=mf_fuel.*burntime;
    oxmass=mf_ox.*burntime;
    lambdas=(fuelmass+oxmass)./(mis+fuelmass+oxmass);
    totalmasses=fuelmass+oxmass+mis;
end