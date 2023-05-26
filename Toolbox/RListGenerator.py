#This should hold all the different ways we generate basic Rlists initially if we ever decide to change it up
#This should be able to use the same parameters as the solidowrks model
import numpy as np
import math
import matplotlib.pyplot as plt

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

def roundRList(xlist, xns, rc, xt, rt, xe, re, rcf, rtf): #rcf = radius_chamber fillet, rtf = radius_throat fillet, 
    
    conslope = ((rc - rt) / (xt - xns))
    divslope = ((re - rt) / (xe - xt)) #slopes are absolute values

    #theta = converging angle, phi = diverging angle

    theta = math.atan(conslope)
    phi = math.atan(divslope)

    #eq1p = rc - rcf 
    #eq1' shows equation of chamber wall shifted down by rcf
    #eq2p = rc - conslope * (x - xns) - (rcf / math.cos(theta)) 
    #eq2' shows equation of converging wall shifted left by rcf*sin(theta) and down by rcf*cos(theta)
    #eq3p = rc - conslope * (x - xns - rtf * math.sin(theta)) + rtf * math.cos(theta) 
    #eq3' shows equation of converging wall shifted right by rtf*sin(theta) and up by rtf*cos(theta)
    #eq4p = rt + divslope * (x - xt + rtf * math.sin(phi)) + rtf * math.cos(phi) 
    #eq4' shows equation of diverging wall shifted left by rtf*sin(phi) and up by rtf*cos(phi)

    #solving for the intersection of two lines is simply m1x+b1=m2x+b2 -> x=(b1-b2)/(m2-m1)
    
    #rcfx and rcfy are the center of the cirlce of the chamber fillet, found at the intersction of eq1 and eq2
    rcfx = ((rc-rcf)-(rc - conslope * ( - xns) - (rcf / math.cos(theta)) ))/(-conslope)
    rcfy = rc-rcf

    rtfx = (( rc - conslope * (- xns ) + rtf / math.cos(theta) ) - \
        (rt + divslope * (- xt ) + rtf / math.cos(phi) ))/\
            (divslope+conslope)
    rtfy = rt + divslope * (rtfx - xt ) + rtf / math.cos(phi) 


    def chamberfillet(rcfx, rcfy, rcf, x): #equation of circle: chamber fillet
        ycf = rcfy + (rcf ** 2 - (x - rcfx) ** 2) ** 0.5 #haven't decided which sign to use
        #change x values to list???
        return ycf

    def throatfillet(rtfx, rtfy, rtf, x): #equation of circle: throat fillet
        ytf = rtfy - (rtf ** 2 - (x - rtfx) ** 2) ** 0.5 #need to choose the sign
        return ytf
    
    realrt = throatfillet(rtfx, rtfy, rtf, rtfx)
    rtdiff = realrt-rt
    xtextended = xt+rtdiff/conslope
    xeextended = xtextended + (re-(rc - conslope * (xtextended - xns)))/divslope

    if xlist[-1]<xeextended:
        dx = np.diff(xlist).min()
        extension = np.arange(xlist[-1]+dx,xeextended+dx, dx)
        xlistnew = np.concatenate((xlist,extension))
    else:
        xlistnew = xlist

    
    # define a new eq4 = re - divslope * (xeextended - x ) + rtf / math.cos(phi) 

    rtfx = (( rc - conslope * (- xns ) + rtf / math.cos(theta) ) - \
        (re - divslope * (xeextended  ) + rtf / math.cos(phi) ))/\
            (divslope+conslope) 
    rtfy = re - divslope * (xeextended - rtfx ) + rtf / math.cos(phi) 

    point1 = [rcfx, rc]
    point2 = [rcfx + rcf * math.sin(theta), rcfy + rcf * math.cos(theta)]
    point3 = [rtfx - rtf * math.sin(theta), rtfy - rtf * math.cos(theta)]
    point4 = [rtfx + rtf * math.sin(phi), rtfy - rtf * math.cos(phi)]

    rlist=np.zeros(xlistnew.size)
    for i in range(xlistnew.size):
        x = xlistnew[i]
        if x < point1[0]: #x before first tangential point (pt 1) on chamber fillet
            rlist[i] = rc
        elif x < point2[0]: #x in  two tangential points (pt 1 and pt 2) of chamber fillet
            rlist[i] = chamberfillet(rcfx, rcfy, rcf, x) #equation of circle: chamber fillet
        elif x < point3[0]: #x before third tangential point (pt 3) on throat fillet
            rlist[i] = rc - conslope * (x - xns) 
        elif x < point4[0]: #x in between two tangential points (pt 3 and pt 4) of throat fillet,
            rlist[i] = throatfillet(rtfx, rtfy, rtf, x) #equation of circle: throat fillet
        else: #x after fourth tangential point (pt 4) on throat fillet
            rlist[i] = re - divslope * ( xeextended - x)

    return rlist, xlistnew

###bezier curve from https://github.com/torresjrjr/Bezier.py by @torresjrjr

"""Bezier, a module for creating Bezier curves.
Version 1.1, from < BezierCurveFunction-v1.ipynb > on 2019-05-02
"""

#import numpy as np

#__all__ = ["Bezier"]
class Bezier():
    def TwoPoints(t, P1, P2):
        """
        Returns a point between P1 and P2, parametised by t.
        INPUTS:
            t     float/int; a parameterisation.
            P1    numpy array; a point.
            P2    numpy array; a point.
        OUTPUTS:
            Q1    numpy array; a point.
        """

        if not isinstance(P1, np.ndarray) or not isinstance(P2, np.ndarray):
            raise TypeError('Points must be an instance of the numpy.ndarray!')
        if not isinstance(t, (int, float)):
            raise TypeError('Parameter t must be an int or float!')

        Q1 = (1 - t) * P1 + t * P2
        return Q1

    def Points(t, points):
        """
        Returns a list of points interpolated by the Bezier process
        INPUTS:
            t            float/int; a parameterisation.
            points       list of numpy arrays; points.
        OUTPUTS:
            newpoints    list of numpy arrays; points.
        """
        newpoints = []
        #print("points =", points, "\n")
        for i1 in range(0, len(points) - 1):
            #print("i1 =", i1)
            #print("points[i1] =", points[i1])

            newpoints += [Bezier.TwoPoints(t, points[i1], points[i1 + 1])]
            #print("newpoints  =", newpoints, "\n")
        return newpoints

    def Point(t, points):
        """
        Returns a point interpolated by the Bezier process
        INPUTS:
            t            float/int; a parameterisation.
            points       list of numpy arrays; points.
        OUTPUTS:
            newpoint     numpy array; a point.
        """
        newpoints = points
        #print("newpoints = ", newpoints)
        while len(newpoints) > 1:
            newpoints = Bezier.Points(t, newpoints)
            #print("newpoints in loop = ", newpoints)

        #print("newpoints = ", newpoints)
        #print("newpoints[0] = ", newpoints[0])
        return newpoints[0]

    def Curve(t_values, points):
        """
        Returns a point interpolated by the Bezier process
        INPUTS:
            t_values     list of floats/ints; a parameterisation.
            points       list of numpy arrays; points.
        OUTPUTS:
            curve        list of numpy arrays; points.
        """

        if not hasattr(t_values, '__iter__'):
            raise TypeError("`t_values` Must be an iterable of integers or floats, of length greater than 0 .")
        if len(t_values) < 1:
            raise TypeError("`t_values` Must be an iterable of integers or floats, of length greater than 0 .")
        if not isinstance(t_values[0], (int, float)):
            raise TypeError("`t_values` Must be an iterable of integers or floats, of length greater than 0 .")

        curve = np.array([[0.0] * len(points[0])])
        for t in t_values:
            #print("curve                  \n", curve)
            #print("Bezier.Point(t, points) \n", Bezier.Point(t, points))

            curve = np.append(curve, [Bezier.Point(t, points)], axis=0)

            #print("curve after            \n", curve, "\n--- --- --- --- --- --- ")
        curve = np.delete(curve, 0, 0)
        #print("curve final            \n", curve, "\n--- --- --- --- --- --- ")
        return curve


#from Bezier import Bezier #
#import numpy as np

#example of a 5-point Bezier curve with parameter t and a numpy array of inital points points1
"""
t_points = np.arange(0, 1, 0.01) #................................. Creates an iterable list from 0 to 1.
points1 = np.array([[0, 0], [0, 8], [5, 10], [9, 7], [4, 3]]) #.... Creates an array of coordinates.
curve1 = Bezier.Curve(t_points, points1) #......................... Returns an array of coordinates.
"""

def paraRlist (xlist, xns, rc, xt, rt_sharp, xe_cone, re_cone, rcf, rtaf, rtef, thetai, thetae, ar): #throat fillet is defined by two sections, 
    #converging "throat approach fillet radius" and diverging "throat expansion fillet radius", 
    #ar = epsilon = nozzle expansion ratio
    conslope = ((rc - rt_sharp) / (xt - xns))
    divslope_c = ((re_cone - rt_sharp) / (xe_cone - xt)) #slopes are absolute values

    #theta = converging angle, phi = diverging angle
    theta = math.atan(conslope)
    phi = math.atan(divslope_c)

    #solving for the intersection of two lines is simply m1x+b1=m2x+b2 -> x=(b1-b2)/(m2-m1)

    #rcfx and rcfy are the center of the cirlce of the chamber fillet, found at the intersction of eq1 and eq2
    rcfx = ((rc-rcf)-(rc - conslope * ( - xns) - (rcf / math.cos(theta)) ))/(-conslope)
    rcfy = rc-rcf
    
    rtfx = xt # center of the circle of the throat approach fillet = (rtfx, rtafy), throat expansion fillet = (rtfx, rtefy)
    rtafy = rt_sharp + rtaf/math.cos(theta)
    rtefy = rtafy - (rtaf - rtef)

    #inflection point I
    ix = rtfx + (rtef * math.sin (thetai)) # center throat expansion fillet + x component 
    iy = rtefy - (rtef * math.sin (thetai))

    lambdaa = 0.5 * (1 + math.cos(phi)) #theoretical correction factor
    ln = (.8 * ((ar**0.5 - 1)*(rt_sharp))) / math.tan (phi) #nozzle length
   
    #exit point E
    exitx = xt + ln
    exity = re_cone #ar**0.5 * (rt_sharp + rtafy - rtaf) #radius of nozzle exit 

    #intersection of I and E
    c1 = iy - math.tan(thetai) * ix
    c2 = exity - math.tan(thetae) * exitx 
    interx = (c2 - c1)/(math.tan(thetai) - math.tan(thetae))
    intery = (math.tan(thetai) * c2 - math.tan(thetae) * c1)/(math.tan(thetai) - math.tan(thetae))

    def chamberfillet(rcfx, rcfy, rcf, x): #equation of circle: chamber fillet
        ycf = rcfy + (rcf ** 2 - (x - rcfx) ** 2) ** 0.5 
        #change x values to list???
        return ycf
    
    ##def throatfillet(rtfx, rtfy, rtf, x): #equation of circle: throat fillet
     #   ytf = rtfy - (rtf ** 2 - (x - rtfx) ** 2) ** 0.5 
     #   return ytf
      
    def throat_a_fillet(rtfx, rtafy, rtaf, x): #equation of circle: throat approach fillet
        ytaf = rtafy - (rtaf ** 2 - (x - rtfx) ** 2) ** 0.5 
        return ytaf

    def throat_e_fillet(rtfx, rtefy, rtef, x): #equation of circle: throat expansion fillet
        ytef = rtefy - (rtef ** 2 - (x - rtfx) ** 2) ** 0.5 
        return ytef

    realrt = throat_a_fillet(rtfx, rtafy, rtaf, rtfx)
    rtdiff = realrt - rt_sharp
    xtextended = xt+rtdiff/conslope
    

    rtfx = xtextended
    rtafy = rt_sharp + rtaf
    rtefy = rtafy - (rtaf - rtef)

    #inflection point I
    ix = rtfx + (rtef * math.sin (thetai)) # center throat expansion fillet + x component 
    iy = rtefy - (rtef * math.cos (thetai))

    lambdaa = 0.5 * (1 + math.cos(phi)) #theoretical correction factor
    ln = (.8 * ((ar**0.5 - 1)*(rt_sharp))) / math.tan (phi) #nozzle length
   
    #exit point E
    exitx = xtextended + ln
    exity = re_cone #ar**0.5 * (rt_sharp + rtafy - rtaf) #radius of nozzle exit 

    #intersection of I and E
    c1 = iy - math.tan(thetai) * ix
    c2 = exity - math.tan(thetae) * exitx 
    interx = (c2 - c1)/(math.tan(thetai) - math.tan(thetae))
    intery = (math.tan(thetai) * c2 - math.tan(thetae) * c1)/(math.tan(thetai) - math.tan(thetae))


    dx = np.diff(xlist).min()
    xlistnew = np.arange(0,exitx,dx)
    point1 = [rcfx, rc]
    point2 = [rcfx + rcf * math.sin(theta), rcfy + rcf * math.cos(theta)]
    point3 = [rtfx - rtaf * math.sin(theta), rtafy - rtaf * math.cos(theta) ] #tangential point to throat approach fillet

    rlist=np.zeros(xlistnew.size)
    t_points = np.linspace(0,1,50000)
    #t_points = np.linspace(ix,exitx,int((exitx-ix)/dx))
    #t_points = np.arange(0, 1, 0.01) #................................. Creates an iterable list from 0 to 1.
    points1 = np.array([[ix, iy], [interx, intery], [exitx, exity]]) #.... Creates an array of coordinates.
    curve1 = Bezier.Curve(t_points, points1) #......................... Returns an array of coordinates.
    
    for i in range(xlistnew.size):
        x = xlistnew[i]
        if x < point1[0]: #x before first tangential point (pt 1) on chamber fillet
            rlist[i] = rc
        elif x < point2[0]: #x in  two tangential points (pt 1 and pt 2) of chamber fillet
            rlist[i] = chamberfillet(rcfx, rcfy, rcf, x) #equation of circle: chamber fillet
        elif x < point3[0]: #x before third tangential point (pt 3) on throat fillet
            rlist[i] = rc - conslope * (x - xns)
        elif x < rtfx: #x from throat approach fillet tangent point to throat
            rlist[i] = throat_a_fillet(rtfx, rtafy, rtaf, x) #y value of throat approach fillet
        elif x < ix: #x from throat to inflection point
            rlist[i] = throat_e_fillet(rtfx, rtefy, rtef, x)#: #y value of throat expansion fillet
        else: #x from inflection point I to exit point E
            tguess = (x-.0001-ix)/(exitx-ix) # prevents x=exitx causing an issue
            rlist[i] = curve1[int(tguess*50000),1]#parabola: Bezier curve
            #rlist[i] = curve1[np.argmin(t_points-x),1]
    #plt.plot([ix,interx,exitx],[iy,intery,exity],'r*')
    return rlist, xlistnew

"""
#the parabolic part was tested and it worked!!
import math 

thetai = 3.14/6
thetae = 8.5*3.14/180

#inflection point I
ix = 0 # center throat expansion fillet + x component 
iy = 4

#exit point E
exitx = 10
exity = 6#radius of nozzle exit 

#intersection of I and E
c1 = iy - math.tan(thetai) * ix
c2 = exity - math.tan(thetae) * exitx 
interx = (c2 - c1)/(math.tan(thetai))
intery = (math.tan(thetai) * c2 - math.tan(thetae) * c1)/(math.tan(thetai) - math.tan(thetae))

t_points = np.arange(0, 1, 0.01) #................................. Creates an iterable list from 0 to 1.
points1 = np.array([[ix, iy], [interx, intery], [exitx, exity]]) #.... Creates an array of coordinates.
curve1 = Bezier.Curve(t_points, points1) #......................... Returns an array of coordinates.

import matplotlib.pyplot as plt
plt.plot(curve1[:,0],curve1[:,1],'r')
plt.plot(points1[:,0],points1[:,1],'go')
plt.show()
"""

