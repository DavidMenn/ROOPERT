"""for determining altitude reached or neccessary rocket to reach target altitude
works by integrating m*dv = (thrust - drag - gravity) dt
Current four function:

rocketEquation(isp, mi = None, thrust = None, burntime = None, L = None, H = None, dt=None, Af=None):
    If you enter H and leave L and mi blank, it returns L to get to H
    if you enter L or mi and leave H blank, returns height reached
    regardles also returns hlist and vlist which define the trajectory

rocketEquationCEA(params, mi = None, thrust = None, burntime = None, L = None, H = None, dt=None, Af=None, ispcorrection = None):
    Same as rocket equation except you pass it params after running first order calcs
    it uses the isp from that to get alittude
    SHOUDL REALLY JUST DO THIS WITH ISENTOPIC EQUATIONS IDIOT

rocketEquationCEA_MassAprox(params, impulse, thrustToWeight, H, dp, dt=None, Af=None, ispcorrection = None):
    Same as above except you can't enter an intert mass or L as that is calcualted from
    the mass approx, so all you do is enter the estimated egnine params and itll give you the
    rocket that makes it there.
    IMPORANT: impulse is an estimate! that might not be enough impulse to get to the alittude,
    so this thing solves for the required impulse to get to that altitude
    thrusttoweight is not an estimate and will be enforced!


ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time= None, title = None ,dt = None):
    plots outputs of rocketequationCEA and CEAMassapprox
"""
import numpy as np
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from ambiance import Atmosphere
import Components.StructuralApproximation as SA
from Toolbox import PressureDropCalculator as PD
heights = np.linspace(-5e3, 80e3, num=1000) #https://pypi.org/project/ambiance/
atmosphere = Atmosphere(heights)

PressureInterpolator = interpolate.interp1d(np.concatenate((heights,[50000000]),axis=0),
    np.concatenate((atmosphere.pressure,[0]),axis=0), kind='linear')
DensityInterpolator = interpolate.interp1d(np.concatenate((heights,[50000000]),axis=0),
    np.concatenate((atmosphere.density,[0]),axis=0), kind='linear')
CdInterpolator = interpolate.interp1d([0,1,2,3,4,5,6,7,8,1000],[.1,.5,.4,.3,.15,.12,.10,.09,.08,0],kind='cubic')

"""Currently two working configs:

Config 1 : h, hlist, vlist = rocketEquation(isp, mi, thrust, burntime, L, H = None, dt=optional, Af=optional)
calculates the altitude reached in meters by a rocket config. Thrust in newtons!!! You can either input inert mass mi
or L, which is mass fraction mp/(mi+mp). h is final altitude reached, and hlist and vlist can be plotted:

plt.plot(np.arange(0,dt*vlist.size,dt),vlist) 

make sure to have dt defined (default to .05)
for acceleration profile (in m/s^2, divide by 9.81 for gs):

alist=np.diff(vlist)/dt
plt.plot(np.arange(0,dt*alist.size,dt),alist) 

Config 2: determines minimum lambda given thrust and isp and height target:
example: (make sure you are in the roopert directory)

import Toolbox.RocketEquation as RE
L, hlist, vlist = RE.rocketEquation(250, thrust =4.48* 2500, burntime = 50, H=100000)

"""
def rocketEquation(isp, mi = None, thrust = None, burntime = None, L = None, H = None, dt=None, Af=None):
    # configure this for a new "none" each time you want to try to calculate some other thing
    if dt is None:
        dt=.05

    if Af is None:
        Af = math.pi *(12*.0254/2)**2 # frontal area, 12 in rocket

    if H is None:

        expectedTime = int((burntime+180)/dt)
        vlist=np.zeros(expectedTime)
        hlist=np.zeros(expectedTime)
        h = 0
        v = 0

        mdot = thrust / (9.81 * isp)
        Mp = mdot * burntime
        if mi is None:
            mi = (1 / L) * Mp - Mp
        M = mi + Mp
        ind = 0

        if thrust < M*9.81:
            raise Exception(f"Thrust too low (Thrust = {thrust}, M = {M})")

        for t in np.arange(0,burntime+dt/2,dt):
            cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2) )/ 2
            M = M - mdot * dt
            h = h + v * dt
            v = v + 9.81 * isp * math.log((M + mdot * dt) / M, math.e) - 9.81 * dt + D / M * dt;
            vlist[ind] = v
            hlist[ind] = h
            ind = ind + 1

        while v > 0:

            h = h + v * dt
            cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
            v = v - 9.81 * dt + D / M * dt
            vlist[ind] = v
            hlist[ind] = h
            ind = ind + 1

        return h, hlist, vlist

    elif L is None and mi is None:
        L=.5
        h=0
        while h < H or h > H + 1000:
            expectedTime = int((burntime + 180) / dt)
            vlist = np.zeros(expectedTime)
            hlist = np.zeros(expectedTime)
            #if h < H:
            #    L = L + .005
            #else:
            #    L = L - .0005

            L = L + .1*(1-h/H)
            print(L)


            expectedTime = int((burntime+180)/dt)
            vlist=np.zeros(expectedTime)
            hlist=np.zeros(expectedTime)
            h = 0
            v = 0

            mdot = thrust / (9.81 * isp)
            Mp = mdot * burntime
            mi = (1 / L) * Mp - Mp
            M = mi + Mp
            ind = 0

            for t in np.arange(0,burntime+dt/2,dt):
                cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
                D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2) )/ 2
                M = M - mdot * dt
                h = h + v * dt
                v = v + 9.81 * isp * math.log((M + mdot * dt) / M, math.e) - 9.81 * dt + D / M * dt
                vlist[ind] = v
                hlist[ind] = h
                ind = ind + 1

            while v > 0:

                h = h + v * dt
                cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
                D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
                v = v - 9.81 * dt + D / M * dt
                vlist[ind] = v
                hlist[ind] = h
                ind = ind + 1

        return L, hlist, vlist

def ShitPlotter(hlist, vlist, thrustlist, isplist, machlist, time= None, title = None ,dt = None):
#MAKE THIS PLOT VELO, Height, and Acceleration in m, m/sm and m/s^2 on the same plot
#Make this plot velo height and accel in mach, km, and g's
#Make this plot density, pressure as fuinction of altitude pls
# make this plot our Cd as function of mach pls
    if time is None:
        time = dt * hlist.size
    if dt is None:
        dt = .05
    if title is None:
        title = "Trajectory Parameters for Unspecified Rocket"
    alist = np.diff(vlist) / dt/9.81
    fig, axs = plt.subplots(2, 3)
    fig.suptitle(title)
    axs[0, 0].plot(np.arange(0, (np.where(hlist==0)[0][1])*dt-dt/2, dt), hlist[0:np.where(hlist==0)[0][1]], 'g')  # row=0, column=0
    axs[1, 0].plot(np.arange(0, (np.where(hlist==0)[0][1]*dt-dt/2), dt), vlist[0:np.where(hlist==0)[0][1]], 'b')  # row=1, column=0
    axs[0, 1].plot(np.arange(0, np.where(hlist==0)[0][1]*dt-dt/2, dt), alist[0:np.where(hlist==0)[0][1]], 'm')  # row=0, column=0
    axs[1, 1].plot(np.arange(0, np.where(hlist==0)[0][1]*dt-dt/2, dt), machlist[0:np.where(hlist==0)[0][1]], 'k')  # row=1, column=0
    axs[0, 2].plot(np.arange(0, time, dt), thrustlist[0:np.size(np.arange(0, time, dt))],
                  'r')  # row=0, column=1
    axs[1, 2].plot(np.arange(0, time, dt), isplist[0:np.size(np.arange(0, time, dt))],
                  'c')  # row=1, column=1
    
    axs[0, 0].set_title('Height')
    axs[1, 0].set_title('Velo (M/s)')
    axs[0, 1].set_title('Acceleration (gs)')
    axs[1, 1].set_title('Mach')
    axs[0, 2].set_title('Thrust (Newtons)')
    axs[1, 2].set_title('Isp (sec)')



def rocketEquationCEA(params, mi = None, thrust = None, burntime = None, L = None, H = None, dt=None, Af=None, ispcorrection = None):
    # configure this for a new "none" each time you want to try to calculate some other thing
    if dt is None:
        dt=.05

    #if ispcorrection is None:
    ispcorrection = 1 #This is legacy, all isp corrections should be done before you get to this point. Just input the right isp in params lol.
    
    if Af is None:
        Af = math.pi *(12*.0254/2)**2 # frontal area, 12 in rocket

    if H is None:

        expectedTime = int((burntime+280)/dt)

        try:
            thrust=params['numengines']*thrust
            expectedTime = int((burntime+180*params['numengines'])/dt)
        except:
            pass
            
        vlist=np.zeros(expectedTime)
        draglist=np.zeros(expectedTime)
        machlist = np.zeros(expectedTime)
        hlist=np.zeros(expectedTime)
        thrustlist = np.zeros(expectedTime)
        isplist = np.zeros(expectedTime)
        h = 0
        v = 0
        
        mdot = thrust / (9.81 * params['isp'])
        Mp = mdot * burntime
        if mi is None:
            mi = (1 / L) * Mp - Mp
        M = mi + Mp
        ind = 0

        if thrust < M*9.81:
            raise Exception(f"Thrust too low (Thrust = {thrust}, M = {M})")

        for t in np.arange(0,burntime+dt/2,dt):
            machlist[ind]=v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h))
            cdtemp = CdInterpolator(machlist[ind])
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2) )/ 2
            M = M - mdot * dt
            h = h + v * dt
            v = v + 9.81 *ispcorrection* params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'], Pamb=PressureInterpolator(h),
                                           frozen=0, frozenAtThroat=0)[0] * math.log((M + mdot * dt) / M, math.e) - 9.81 * dt + D / M * dt
            thrustlist[ind] = 9.81 * ispcorrection* params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'], Pamb=PressureInterpolator(h),
                                           frozen=0, frozenAtThroat=0)[0] * mdot
            isplist[ind] = ispcorrection*params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'], Pamb=PressureInterpolator(h),
                                           frozen=0, frozenAtThroat=0)[0]
            vlist[ind] = v
            hlist[ind] = h
            draglist[ind] = D
            ind = ind + 1
            if ind%100 == 0:
                print(t)
        while v > 0:
            
            h = h + v * dt
            machlist[ind] = v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h))
            cdtemp = CdInterpolator(machlist[ind])
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
            v = v - 9.81 * dt + D / M * dt
            vlist[ind] = v
            hlist[ind] = h
            draglist[ind] = D
            ind = ind + 1
        print(f"max height reached = {h} meters")
        return h, hlist, vlist, thrustlist, isplist, machlist

    elif L is None and mi is None:
        L=.5
        h=0
        while h < H or h > H + 1000:
            expectedTime = int((burntime + 180) / dt)
            vlist = np.zeros(expectedTime)
            hlist = np.zeros(expectedTime)
            machlist = np.zeros(expectedTime)
            thrustlist = np.zeros(expectedTime)
            isplist = np.zeros(expectedTime )
            #if h < H:
            #    L = L + .005
            #else:
            #    L = L - .0005

            L = L + .1*(1-h/H)
            print(L)


            expectedTime = int((burntime+180)/dt)
            vlist=np.zeros(expectedTime)
            hlist=np.zeros(expectedTime)
            h = 0
            v = 0

            mdot = thrust / (9.81 * params['isp'])
            Mp = mdot * burntime
            mi = (1 / L) * Mp - Mp
            M = mi + Mp
            ind = 0
            print(f"M = {M}")
            if thrust < M * 9.81:
                raise Exception(f"Thrust too low (Thrust = {thrust}, M = {M})")
            for t in np.arange(0,burntime+dt/2,dt):
                machlist[ind] = v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h))
                cdtemp = CdInterpolator(machlist[ind])
                D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2) )/ 2
                M = M - mdot * dt
                h = h + v * dt
                v = v + 9.81 *ispcorrection* params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'],
                                                                  Pamb=PressureInterpolator(h),
                                                                  frozen=0, frozenAtThroat=0)[0] * math.log(
                    (M + mdot * dt) / M, math.e) - 9.81 * dt + D / M * dt
                thrustlist[ind] = 9.81 * \
                                  ispcorrection*params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'],
                                                                     Pamb=PressureInterpolator(h),
                                                                     frozen=0, frozenAtThroat=0)[0] * mdot
                isplist[ind] = ispcorrection*params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'],
                                                                  Pamb=PressureInterpolator(h),
                                                                  frozen=0, frozenAtThroat=0)[0]
                vlist[ind] = v
                hlist[ind] = h
                ind = ind + 1

            while v > 0:

                h = h + v * dt
                machlist[ind] = v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h))
                cdtemp = CdInterpolator(machlist[ind])
                D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
                v = v - 9.81 * dt + D / M * dt
                vlist[ind] = v
                hlist[ind] = h
                ind = ind + 1
        print(f"max height reached = {h} meters")
        return L, hlist, vlist, thrustlist, isplist, machlist

def rocketEquationCEA_MassAprox(params, impulse, thrustToWeight, H, dp, dt=None, Af=None, ispcorrection = None):
    # configure this for a new "none" each time you want to try to calculate some other thing
    hslist=np.zeros(1000)
    impulselist=np.zeros(1000)
    if dt is None:
        dt=.05

    if ispcorrection is None:
        ispcorrection = 1

    if Af is None:
        Af = math.pi *(12*.0254/2)**2 # frontal area, 12 in rocket
    miinittemp = 0
    mi=100000
    totalmass = 0
    #You gave me an impulse estimate and required thrust to weight,
    # so im finding the rocket with that impulse and thurst to weight
    while abs(miinittemp-mi)>2:
        miinittemp=mi
        mi, L, totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank = SA.mass_approx(params['pc'], params['dp'], 12, params['rho_fuel'], params['rho_ox'],
                                                       params['thrust'], params['isp'], params['time'], params['rm'])
        params['thrust'] = totalmass * thrustToWeight * 9.81  # recompute with a resonable thrust
        params['time'] = impulse/params['thrust']
        params['mdot'] = params['thrust']/params['isp']/params['g']
        params['mdot_fuel'] = params['mdot']/(params['rm']+1)
        params['mdot_ox'] = params['rm']*params['mdot']/(params['rm']+1)
        oxpd, fuelpd, dparray = PD.pressureDrop(OxDensityMetric = params['rho_ox'],     # kg/m^3
                                    OxMassFlowMetric = params['mdot_ox'], # kg/s
                                    OxTubeDiam = 3/4, #in  
                                    FuelDensityMetric = params['rho_fuel'],       # kg/m^3
                                    FuelMassFlowMetric = params['mdot_fuel'],# kg/s
                                    FuelTubeDiam = 3/4 ,
                                    params = params               )
        params['dp'] = dparray 
        print(f"MI is {mi} and Thrust is {params['thrust']}")
    h=H-1
    impulse1 = impulse+10000
    impulse2 = impulse+1000
    h1=h+100
    h2=h+10
    index=-1
    while h < H+10 or h > H + 500:
        index=index+1
        #if h < H:
        #    L = L + .005
        #else:
        #    L = L - .0005
        print(f"got to {h} at impulse {impulse} and mi = {mi}")
        #impulse = impulse*(H/h)
        #dhdimp=(h2-h)/(impulse2-impulse)
        #dhdimp2=((h1-h-(impulse1-impulse)*dhdimp)/((impulse1-impulse)**2)*2+ (h2-h-(impulse2-impulse)*dhdimp)/((impulse2-impulse)**2)*2)/2
        #impulse1 = impulse2
        #h1 = h2
        #impulse2 = impulse
        #h2=h
        #if h<H:

            #deltaimpulse = np.max(np.roots([dhdimp2/2, dhdimp, h-H]))
        #else:

            #deltaimpulse = np.max(np.roots([dhdimp2 / 2, dhdimp, h - H]))
        #deltaimpulse = np.max(np.roots([ dhdimp, h-H]))/2.5 + ((np.max(np.roots([ dhdimp, h-H]))>0)*2-1)*50
        deltaimpulse = (impulse*H/h-impulse)*.75
        impulse = impulse+deltaimpulse
        miinittemp = 0
        while abs(miinittemp - mi) > 2:
            miinittemp = mi
            params['thrust'] = totalmass * thrustToWeight * 9.81  # recompute with a resonable thrust
            params['time'] = impulse/params['thrust']
            params['mdot'] = params['thrust']/params['isp']/params['g']
            params['mdot_fuel'] = params['mdot']/(params['rm']+1)
            params['mdot_ox'] = params['rm']*params['mdot']/(params['rm']+1)
            oxpd, fuelpd, dparray = PD.pressureDrop(OxDensityMetric = params['rho_ox'],     # kg/m^3
                                        OxMassFlowMetric = params['mdot_ox'], # kg/s
                                        OxTubeDiam = 3/4, #in  
                                        FuelDensityMetric = params['rho_fuel'],       # kg/m^3
                                        FuelMassFlowMetric = params['mdot_fuel'],# kg/s
                                        FuelTubeDiam = 3/4 ,
                                        params = params               )
            params['dp'] = dparray 
            mi, L, totalmass, wstruct, newheight, heightox, heightfuel, vol_ox, vol_fuel, P_tank = SA.mass_approx(params['pc'], params['dp'], 12, params['rho_fuel'], params['rho_ox'],
                                              params['thrust'], params['isp'], params['time'], params['rm'])

            print(f"MI is {mi}, L is {L} and Thrust is {params['thrust']} and dp is {params['dp']}")

        expectedTime = int((params['time'] + 180) / dt)+1000
        vlist = np.zeros(expectedTime)
        draglist = np.zeros(expectedTime)
        hlist = np.zeros(expectedTime)
        machlist = np.zeros(expectedTime)
        thrustlist = np.zeros(expectedTime)
        isplist = np.zeros(expectedTime)
        h = 0
        v = 0

        mdot = params['thrust'] / (9.81 * params['isp'])
        Mp = mdot * params['time']
        M = mi + Mp
        ind = 0
        print(f"M = {M}")
        if params['thrust'] < M * 9.81:
            raise Exception(f"Thrust too low (Thrust = {params['thrust']}, M = {M})")
        for t in np.arange(0,params['time']+dt/2,dt):
            machlist[ind] = v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h))
            cdtemp = CdInterpolator(machlist[ind])
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2) )/ 2
            M = M - mdot * dt
            h = h + v * dt
            ispamb = params['CEA'].estimate_Ambient_Isp(Pc=params['pc'], MR=params['rm'], eps=params['er'],
                                                              Pamb=PressureInterpolator(h),
                                                              frozen=0, frozenAtThroat=0)[0]
            v = v + 9.81 *ispcorrection* ispamb * math.log(
                (M + mdot * dt) / M, math.e) - 9.81 * dt + D / M * dt
            thrustlist[ind] = 9.81 * \
                              ispcorrection*ispamb* mdot
            isplist[ind] = ispcorrection*ispamb
            vlist[ind] = v
            hlist[ind] = h
            draglist[ind] = D
            ind = ind + 1

        while v > 0:

            h = h + v * dt
            machlist[ind] = v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h))
            cdtemp = CdInterpolator(machlist[ind])
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
            v = v - 9.81 * dt + D / M * dt
            vlist[ind] = v
            hlist[ind] = h
            draglist[ind] = D
            ind = ind + 1
        impulselist[index]=impulse
        hslist[index]=h
    print(f"got to {h} at impulse {impulse} and mi = {mi}")
    print(f"MI is {mi}, L is {L} and Thrust is {params['thrust']}")
    #plt.plot(np.arange(0,60,dt),draglist[0:int(60/dt)], label="Drag")
    ##plt.plot(np.arange(0,60,dt),thrustlist[0:int(60/dt)], label="Thrust")
    #plt.ylabel("Force [N]")
    #plt.xlabel("Time [s]")
    #plt.legend()
    #plt.title("Forces on rocket for first minute of flight")
    return L, mi, hlist, vlist, thrustlist, isplist, machlist#, hslist, impulselist
