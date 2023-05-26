#For eventual flow calcs to make them more legible
import numpy as np
import math
import scipy
def totalT(T, gam, M):
    return T * (1 + ((gam - 1) / 2) * (M ** 2))

def TFromTotalT(totalT, gam, M):
    return totalT / (1 + ((gam - 1) / 2) * (M ** 2))

def totalP(P, gam, M):
    return P * (1 + ((gam - 1) / 2) * (M ** 2))**(gam/(gam-1))

def PFromTotalP(totalP, gam, M):
    return totalP * (1 + ((gam - 1) / 2) * (M ** 2))**(-gam/(gam-1))

def machFromP(totalP,P,gam):
    return (((totalP/P)**((gam-1)/gam)-1)/((gam - 1) / 2))**.5
#Going to use this in the future to figure out injector tolerances
def PFromMdotAtThroat(mdot, gam, totalT, at, R):
    totalP = mdot * (totalT ** .5) / at / ((gam / R) ** .5) * (((gam + 1) / 2) ** ((gam + 1) / (2 * (gam - 1))))
    return PFromTotalP(totalP, gam, 1)



def AreaForChokedFlow(pt,Tt,gam,mdot,specificR):

    return mdot * math.sqrt(Tt) / pt / math.sqrt(gam / specificR) / (
                (1 + (gam - 1) / 2) ** (-(gam + 1) / (2 * (gam - 1))))

def RtoK(degR):
    return degR*5/9

def machFromArea(a,at,gam,supersonic = False):
    machtol=.0001
    gm1 = gam - 1
    gp1 = gam + 1
    ARatio = a/at
 
    #eq = lambda M : abs(at/a - M * ((gp1/2) / (1+(gm1/2)* M**2))**(gp1/gm1/2))  
    eq = lambda M : abs((1/M**2)*(((2+gm1*M**2)/gp1)**(gp1/gm1))-ARatio**2)
 
    if not supersonic:
        return scipy.optimize.minimize_scalar(eq, bounds=(1e-6, 1), tol=machtol, method='bounded')['x']
    else:
        return scipy.optimize.minimize_scalar(eq, bounds=(1+1e-6, 50), tol=machtol, method='bounded')['x']


    """
    if supersonic:
        P = 2 / (gam + 1)
        Q = 1 - P
        E = 1 / Q
        R = (a / at) ** (2*Q/P)
        A = Q ** (1 / P)
        r = (R - 1) / (2 * A)
        X = 1 / (1 + r + (r * (r + 2)) ** .5)
        Xnew = P * (X - 1) / (1 - R * (P + Q * X) ** (-P / Q))
        tol = ((1/machtol)**2)
        while abs(X-Xnew)>tol:
            X=Xnew
            Xnew = P * (X - 1) / (1 - R * (P + Q * X) ** (-P / Q))
        mach = 1/(Xnew**.5)
    else:
        P = 2 / (gam + 1)
        Q = 1 - P
        E = 1 / Q
        R = (a / at) ** (2)
        A = P ** (1 / Q)
        r = (R - 1) / (2 * A)
        X = 1 / (1 + r + (r * (r + 2)) ** .5)
        Xnew = P * (X - 1) / (1 - R * (P + Q * X) ** (-P / Q))
        tol = ((machtol) ** 2)
        while abs(X - Xnew) > tol:
            X = Xnew
            Xnew = P * (X - 1) / (1 - R * (P + Q * X) ** (-P / Q))
        mach = (Xnew ** .5)

    return mach
    """
