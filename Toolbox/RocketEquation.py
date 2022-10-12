#USE ROCKET CEA OR EVEN JUST ISENTROPIC EQUATIONS TO GET ISP AS WE ASCEND YOU DUNCE

import numpy as np
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from ambiance import Atmosphere

#heights = np.linspace(-5e3, 80e3, num=1000) #https://pypi.org/project/ambiance/
#atmosphere = Atmosphere(heights)

PressureInterpolator = interpolate.interp1d([0 ,11000 ,20000 ,32000 ,47000, 51000, 71000, 100000000],
    [101325,22732,5474,868,110,66,3, 0], kind='linear')
DensityInterpolator = interpolate.interp1d([0 ,11000 ,20000 ,32000 ,47000, 51000, 71000, 100000000],
    [1.225,.36391,.08803,.01322,.00143,.00086,.000064,0], kind='linear')
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

def ShitPlotter(hlist, vlist, dt = None):
#MAKE THIS PLOT VELO, Height, and Acceleration in m, m/sm and m/s^2 on the same plot
#Make this plot velo height and accel in mach, km, and g's
#Make this plot density, pressure as fuinction of altitude pls
# make this plot our Cd as function of mach pls
    if dt is None:
        dt = .05
    alist = np.diff(vlist) / dt
    plt.plot(np.arange(0, dt * alist.size, dt), alist, )
    fig, axs = plt.subplots(2, 2)
    fig.suptitle('Sharing x per column, y per row')
    axs[0, 0].plot(x, y)
    axs[0, 0].set_title('Axis [0, 0]')
    axs[0, 1].plot(x, y, 'tab:orange')
    axs[0, 1].set_title('Axis [0, 1]')
    axs[1, 0].plot(x, -y, 'tab:green')
    axs[1, 0].set_title('Axis [1, 0]')
    axs[1, 1].plot(x, -y, 'tab:red')
    axs[1, 1].set_title('Axis [1, 1]')

    for ax in axs.flat:
        ax.set(xlabel='x-label', ylabel='y-label')

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
            cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
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
            cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
            D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
            v = v - 9.81 * dt + D / M * dt
            vlist[ind] = v
            hlist[ind] = h
            ind = ind + 1

        return h, hlist, vlist, thrustlist, isplist

    elif L is None and mi is None:
        L=.5
        h=0
        while h < H or h > H + 1000:
            expectedTime = int((burntime + 180) / dt)
            vlist = np.zeros(expectedTime)
            hlist = np.zeros(expectedTime)
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

            for t in np.arange(0,burntime+dt/2,dt):
                cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
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
                cdtemp = CdInterpolator(v / math.sqrt(1.4 * PressureInterpolator(h) / DensityInterpolator(h)))
                D = -(cdtemp * Af * DensityInterpolator(h) * (v ** 2)) / 2
                v = v - 9.81 * dt + D / M * dt
                vlist[ind] = v
                hlist[ind] = h
                ind = ind + 1

        return L, hlist, vlist, thrustlist, isplist
