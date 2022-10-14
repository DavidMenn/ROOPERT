#USE ROCKET CEA OR EVEN JUST ISENTROPIC EQUATIONS TO GET ISP AS WE ASCEND YOU DUNCE

import numpy as np
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from ambiance import Atmosphere
import Components.StructuralApproximation as SA
heights = np.linspace(-5e3, 80e3, num=1000) #https://pypi.org/project/ambiance/
atmosphere = Atmosphere(heights)

PressureInterpolator = interpolate.interp1d(np.concatenate((heights,[500000]),axis=0),
    np.concatenate((atmosphere.pressure,[0]),axis=0), kind='linear')
DensityInterpolator = interpolate.interp1d(np.concatenate((heights,[500000]),axis=0),
    np.concatenate((atmosphere.density,[0]),axis=0), kind='linear')
CdInterpolator = interpolate.interp1d([0,1,2,3,4,5,6,7,8],[.1,.5,.4,.3,.15,.12,.10,.09,.08],kind='cubic')

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

    axs[0, 0].plot(np.arange(0, dt * hlist.size, dt), hlist, 'g')  # row=0, column=0
    axs[1, 0].plot(np.arange(0, dt * hlist.size, dt), vlist, 'b')  # row=1, column=0
    axs[0, 1].plot(np.arange(0, dt * alist.size, dt), alist, 'g')  # row=0, column=0
    axs[1, 1].plot(np.arange(0, dt * hlist.size, dt), machlist, 'b')  # row=1, column=0
    axs[0, 2].plot(np.arange(0, time, dt), thrustlist[0:np.size(np.arange(0, time, dt))],
                  'r')  # row=0, column=1
    axs[1, 2].plot(np.arange(0, time, dt), isplist[0:np.size(np.arange(0, time, dt))],
                  'k')  # row=1, column=1
    
    axs[0, 0].set_title('Height')
    axs[1, 0].set_title('Velo (M/s)')
    axs[0, 1].set_title('Acceleration (gs)')
    axs[1, 1].set_title('Mach')
    axs[0, 2].set_title('Thrust (Newtons)')
    axs[1, 2].set_title('Isp (sec)')
    plt.show()

    for ax in axs.flat:
        ax.set(xlabel='time')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

def rocketEquationCEA(params, mi = None, thrust = None, burntime = None, L = None, H = None, dt=None, Af=None, ispcorrection = None):
    # configure this for a new "none" each time you want to try to calculate some other thing
    if dt is None:
        dt=.05

    if ispcorrection is None:
        ispcorrection = 1

    if Af is None:
        Af = math.pi *(12*.0254/2)**2 # frontal area, 12 in rocket

    if H is None:

        expectedTime = int((burntime+180)/dt)
        vlist=np.zeros(expectedTime)
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

        return L, hlist, vlist, thrustlist, isplist, machlist

def rocketEquationCEA_MassAprox(params, impulse, thrustToWeight, H, dp, dt=None, Af=None, ispcorrection = None):
    # configure this for a new "none" each time you want to try to calculate some other thing
    if dt is None:
        dt=.05

    if ispcorrection is None:
        ispcorrection = 1

    if Af is None:
        Af = math.pi *(12*.0254/2)**2 # frontal area, 12 in rocket
    miinittemp = 0
    mi=100000
    totalmass = 0
    while abs(miinittemp-mi)>2:
        miinittemp=mi
        mi, L, totalmass = SA.mass_approx(params['pc'], dp, 12, params['rho_fuel'], params['rho_ox'],
                                                       params['thrust'], params['isp'], params['time'], params['rm'])
        params['thrust'] = totalmass * thrustToWeight * 9.81  # recompute with a resonable thrust
        params['time'] = impulse/params['thrust']
        print(f"MI is {mi} and Thrust is {params['thrust']}")
    h=H-1
    while h < H or h > H + 500:

        #if h < H:
        #    L = L + .005
        #else:
        #    L = L - .0005
        print(f"got to {h} at impulse {impulse} and mi = {mi}")
        impulse = impulse*(H/h)
        miinittemp = 0
        while abs(miinittemp - mi) > 2:
            miinittemp = mi
            params['thrust'] = totalmass * thrustToWeight * 9.81  # recompute with a resonable thrust
            params['time'] = impulse / params['thrust']
            mi, L, totalmass = SA.mass_approx(params['pc'], dp, 12, params['rho_fuel'], params['rho_ox'],
                                              params['thrust'], params['isp'], params['time'], params['rm'])

            print(f"MI is {mi}, L is {L} and Thrust is {params['thrust']}")

        expectedTime = int((params['time'] + 180) / dt)
        vlist = np.zeros(expectedTime)
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

    return L, mi, hlist, vlist, thrustlist, isplist, machlist
