#This should hold all the different ways we generate basic Rlists initially if we ever decide to change it up
import numpy as np


def sharpRList(xlist, xns, rc, xt, rt, xe, re):
    rlist=np.zeros(xlist.size)
    conslope = ((rc - rt) / (xt - xns))
    divslope = ((re - rt) / (xe - xt))
    for i in range(xlist.size):
        x = xlist[i]
        if x < xns:
            rlist[i] = rc
        elif x < xt:
            rlist[i] = rc - conslope * (x - xns)
        elif x == xt:
            rlist[i] = rt
        else:
            rlist[i] = rt + divslope * (x - xt)

    return rlist