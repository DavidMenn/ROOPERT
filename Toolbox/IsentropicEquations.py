import numpy as np
import math
def totalT(T, gam, M):
    return T * (1 + ((gam - 1) / 2) * (M ** 2))

def TFromTotalT(totalT, gam, M):
    return totalT / (1 + ((gam - 1) / 2) * (M ** 2))

def totalP(P, gam, M):
    return P * (1 + ((gam - 1) / 2) * (M ** 2))**(gam/(gam-1))

def PFromTotalP(totalP, gam, M):
    return totalP * (1 + ((gam - 1) / 2) * (M ** 2))**(-gam/(gam-1))

def PFromMdotAtThroat(mdot, gam, totalT, at, R):
    totalP = mdot * (totalT ** .5) / at / ((gam / R) ** .5) * (((gam + 1) / 2) ** ((gam + 1) / (2 * (gam - 1))))
    return PFromTotalP(totalP, gam, 1)



def AreaForChokedFlow(p,T,gam,mdot,specificR):
    Tt=totalT(T,gam,1)
    pt=totalP(p,gam,1)
    return mdot * math.sqrt(Tt) / pt / math.sqrt(gam / specificR) / (
                (1 + (gam - 1) / 2) ** (-(gam + 1) / (2 * (gam - 1))))

def RtoK(degR):
    return degR*5/9
