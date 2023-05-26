import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import optimize
def beanfunc(x,a,b,c,d,s, sign=1): # here x is in the theta direction (basiswidth), y is in the r direction (basisheight)
    # bean func is x^2 + (y-d)^2 = (a(y-d)^2+bx^2+c(y-d))
    homogenousbean = lambda y : (a*(y*s-d)**2+b*x**2+c*(y*s-d))**2 -x**2 - (y*s-d)**2
    result = optimize.root(homogenousbean, sign*10, args=(), method='hybr', jac=None, tol=None, callback=None, options=None)
    return result
def ellipsefunc(x,axisy,axisx = .5, sign=1):
    return sign*(1-(x/axisx)**2)**.5*axisy
def area(beanpoints): #beanpoints of the form [0,.1,...,1 ; topleft, ... , topright; bottomleft, ..., bottomright]
    sum = 0
    for i in range(1,beanpoints.shape[1]):
        dx = abs(beanpoints[0,i]-beanpoints[0,i-1])
        sum = sum + abs(beanpoints[1,i]*dx) + abs(beanpoints[2,i]*dx)
    return sum
def perimeter(beanpoints):
    sum = beanpoints[1,0]+beanpoints[1,-1]+beanpoints[2,0]+beanpoints[2,-1]
    for i in range(1,beanpoints.shape[1]):
        dx = abs(beanpoints[0,i]-beanpoints[0,i-1])
        sum = sum + np.linalg.norm([beanpoints[1,i]-beanpoints[1,i-1],dx]) + np.linalg.norm([beanpoints[2,i]-beanpoints[2,i-1],dx])
    return sum
def beanmaker(axisyval,ch,beanpoints,a,b,c,d,s,areadesired,cw): 
    for index in range(0,beanpoints.shape[1]):
        beanpoints[1,index] = ellipsefunc(beanpoints[0,index],axisy=axisyval,sign=1)*ch
        beanpoints[2,index] = beanfunc(beanpoints[0,index],a,b,c,d,s,sign=-1).x[0]*ch
    beanpoints[0,:] = beanpoints[0,:]*cw
    #plt.plot([cw*.5, cw*.5, -cw*.5, -cw*.5],[areadesired/cw*.5,-areadesired/cw*.5,-areadesired/cw*.5,areadesired/cw*.5,])
    #plt.plot(beanpoints[0,:],beanpoints[1,:])
    ##plt.plot(beanpoints[0,:],beanpoints[2,:])
    #plt.show()
    #return np.linalg.norm([abs(area(beanpoints)-areadesired)/areadesired,abs(perimeter(beanpoints)-perimdesired)/perimdesired])
    return area(beanpoints),perimeter(beanpoints),beanpoints

numbeanpoints = 40
beanpoints = np.linspace(-.5,.5,numbeanpoints)
beanpoints = np.vstack((beanpoints,np.zeros((2,numbeanpoints))))
beanpoints[1,:] = ellipsefunc(beanpoints[0,:],axisy = .35,axisx = .5, sign=1)+.5
beanpoints[2,:] =  ellipsefunc(beanpoints[0,:],axisy = .4,axisx = .5, sign=-1)+.5 + .175*(np.cos(beanpoints[0,:]*math.pi*2)+1)
beanreduce = np.min(beanpoints[2,:])
beanpoints[1,:] = beanpoints[1,:] - beanreduce
beanpoints[2,:]= beanpoints[2,:] - beanreduce
beanmultiply = np.max(beanpoints[1,:])
beanpoints[1,:] = beanpoints[1,:]/beanmultiply
beanpoints[2,:]= beanpoints[2,:] /beanmultiply
plt.plot(beanpoints[0,:],beanpoints[1,:])
plt.plot(beanpoints[0,:],beanpoints[2,:])
plt.show()

"""areadesired = chlist[i]*cwlist[i]
perimdesired = 2*chlist[i]+2*cwlist[i]
axisyval = .2
chguess = chlist[i]
areacurrent, perim, beanpoints = beanmaker(axisyval,chguess,beanpoints,33,2,.3,.005,.127,areadesired,cwlist[i])
tol=.0001
while abs(areacurrent-areadesired)/areadesired>tol:
    beanpoints = np.linspace(-.5,.5,numbeanpoints)
    beanpoints = np.vstack((beanpoints,np.zeros((2,numbeanpoints))))
    
    chguess = chguess- (areacurrent-areadesired)/areadesired*chguess
    areacurrent, perim, beanpoints = beanmaker(axisyval,chguess,beanpoints,33,2,.3,.005,.127,areadesired,cwlist[i])

#plt.plot(beanpoints[0,:],beanpoints[1,:])
#plt.plot(beanpoints[0,:],beanpoints[2,:])
beanpoints[1,:] = -np.min(beanpoints[2,:])+beanpoints[1,:] #currently bottom is at -.5*ch, gotta shift it up
beanpoints[2,:] = -np.min(beanpoints[2,:]) + beanpoints[2,:]
#plt.plot(beanpoints[0,:],beanpoints[1,:])
#plt.plot(beanpoints[0,:],beanpoints[2,:])
#plt.show()
hydraulicdiamlistnew[i] = 4*areacurrent/perim
chlistnew[i] = np.max(beanpoints[1,:])"""