import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import optimize
def ChanelBoxCorners(xlistoriginal,rlist,twlist,helicitylist, chlist, cwlist):
    ylist = np.zeros(len(xlistoriginal))
    xlist=np.zeros(len(ylist))
    zlist=np.zeros(len(ylist))
    xlist[0]=rlist[0]+twlist[0]
    theta = 0
    for i in range(1,len(ylist)):
        if i%35 == 0: #this is just for debugging
            i=i
        if i == 0:
            drdy = (rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i])
        elif i == len(xlistoriginal)-1:
            drdy = (rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])
        else:
            drdy = ((rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])+(rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i]))/2
        axialangle = -math.atan(drdy)
        ylist[i] = xlistoriginal[i]
        r=rlist[i]+twlist[i]*math.cos(axialangle)
        #c=math.tan(helicitylist[i])*2*math.pi*r
        dy=(ylist[i]-ylist[i-1])
        if helicitylist[i]==math.pi/2:
            dcircum = 0
        else:
            dcircum = dy/math.tan(helicitylist[i])
        dtheta = dcircum/r
        theta=theta+dtheta
        xlist[i]=r*math.cos(theta)
        zlist[i]=r*math.sin(theta)
    # Now we hzve the geometry for the sweep curve, or the curve along rlist+twlist
    #We want to make four lists now, one for each corner of the box
    ylistnew = np.zeros((4,len(xlistoriginal)))
    xlistnew=np.zeros((4,len(ylist)))
    zlistnew =np.zeros((4,len(ylist)))
    #fig = plt.figure() # for showing vectors for debugging
    #ax=fig.add_subplot(projection='3d')
    for i in range(0,len(ylist)):
        angle = abs(math.atan(xlist[i]/zlist[i]))
        #xlist[i]=xlist[i]+math.sin(angle)*np.sign(xlist[i])*chlist[i]
        #zlist[i]=zlist[i]+math.cos(angle)*np.sign(zlist[i])*chlist[i]
        #Find the angle of rlist (converging angle, divering angle, etc)
        if i%35 == 0: #this is just for debugging
            i=i
        if i == 0:
            drdy = (rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i])
        elif i == len(xlistoriginal)-1:
            drdy = (rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])
        else:
            drdy = ((rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])+(rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i]))/2
        axialangle = -math.atan(drdy)
        # Now find the basis vectors for the plane with the curve as its normal
        #basis1 = [math.sin(angle)*np.sign(xlist[i]), 0 , math.cos(angle)*np.sign(zlist[i])]
        #basis1*basis2 = 0
        #basis2*[0 1 0] = cos(90-rlistangle)*|basis2|
        #basis2 = [A*math.cos(angle)*np.sign(zlist[i]), cos(math.pi/2-rlistangle)*|basis2|, -A*math.sin(angle)*np.sign(xlist[i])]
        #Set |basis2| =1
        #sqrt(A^2+cos(math.pi/2-rlistangle)^2)=1
        #A^2=1-cos(math.pi/2-rlistangle)^2
        Aaxial=math.sqrt(1-(math.cos(math.pi/2-axialangle)**2))
        #Ahelicity=math.sqrt(1-(math.cos(helicitylist[i])**2))
        basisheight = [Aaxial*math.sin(angle)*np.sign(xlist[i]), math.cos(math.pi/2-axialangle), Aaxial*math.cos(angle)*np.sign(zlist[i])]
        #now with basis height, we want to find a perpendicular vector that is also angled with the helicity
        # for a helicity of 90, the axial (y) component should be zero
        #chaenlcurve vector defined in cylydnircal coordinates : (r,theta,y)
        #-theta/sin(hel)=y/cos(hel) 
        #r/sin(axial)=y/cos(axial)
        #r^2+theta^2+y^2=1
        #ytan(hel)^2+ytan(axial)^2+y^2=1
        #y=math.sqrt(1/(1+math.tan(math.pi/2-helicitylist[i])**2+math.tan(axialangle)**2))
        # x = r*sin(angle) +theta*cos(angle)
        # z = r*cos(angle) - theta*sin(agnle)
        # chanelvector = [math.tan(axialangle)*y*math.sin(angle) - math.tan(helicitylist[i])*y*math.cos(angle)
        # , y, 
        # math.tan(axialangle)*y*math.cos(angle) + math.tan(helicitylist[i])*y*math.sin(angle)]
        #cylyndrilca basis height(r,theta,z) = (cos(axialangle),0,sin(axialangle))
        #cylidrilca chanel vector = [-tan(axialangle)y, -tan(hel)y, y]
        angle = np.arctan2(xlist[i],zlist[i]) # this is so ham but I did it janky before, now I fix
        if helicitylist[i] == math.pi/2:
            y=math.sqrt(1/(1+math.tan(axialangle)**2))
            chanelvector = [-(math.tan(axialangle)*y*math.sin(angle)),
            y, 
            -(math.tan(axialangle)*y*math.cos(angle))]
        else:
            y=math.sqrt(1/(1+(1/math.tan(helicitylist[i]))**2+math.tan(axialangle)**2))
            chanelvector = [-(math.tan(axialangle)*y*math.sin(angle) +(1/math.tan(helicitylist[i]))*y*math.cos(angle)),
            y, 
            -(math.tan(axialangle)*y*math.cos(angle) - (1/math.tan(helicitylist[i]))*y*math.sin(angle))]
        
        

        # the perpendicular vector is just chanelvector cross basisheight
        basiswidth = np.cross(chanelvector,basisheight)
        if i == 0: # this is for the exhaust-side extension (to ensure complete cutting)
            chanelvectorinit = chanelvector
            basisheightinit = basisheight
            basiswidthinit = basiswidth
        #Helps verify vectors point where they should
        #ax.plot([xlist[i],xlist[i] + basiswidth[0]*.01],[ylist[i],ylist[i] + basiswidth[1]*.01],
        #                [zlist[i],zlist[i] + basiswidth[2]*.01],'r')
        #ax.plot([xlist[i],xlist[i] + chanelvector[0]*.01],[ylist[i],ylist[i] + chanelvector[1]*.01],
        #                [zlist[i],zlist[i] + chanelvector[2]*.01],'g')
        #ax.plot([xlist[i],xlist[i] + basisheight[0]*.01],[ylist[i],ylist[i] + basisheight[1]*.01],
        #                [zlist[i],zlist[i] + basisheight[2]*.01],'b')

        #basiswidth = [Ahelicity*math.cos(angle)*np.sign(zlist[i]), math.cos(helicitylist[i]), -Ahelicity*math.sin(angle)*np.sign(xlist[i])]
        xlistnew[0,i] = xlist[i] + basiswidth[0]*cwlist[i]*.5
        ylistnew[0,i] = ylist[i] + basiswidth[1]*cwlist[i]*.5
        zlistnew[0,i] = zlist[i] + basiswidth[2]*cwlist[i]*.5

        xlistnew[1,i] = xlist[i] + basiswidth[0]*cwlist[i]*.5 + basisheight[0]*chlist[i]
        ylistnew[1,i] = ylist[i] + basiswidth[1]*cwlist[i]*.5 + basisheight[1]*chlist[i]
        zlistnew[1,i] = zlist[i] + basiswidth[2]*cwlist[i]*.5 + basisheight[2]*chlist[i]

        xlistnew[2,i] = xlist[i] - basiswidth[0]*cwlist[i]*.5 + basisheight[0]*chlist[i]
        ylistnew[2,i] = ylist[i] - basiswidth[1]*cwlist[i]*.5 + basisheight[1]*chlist[i]
        zlistnew[2,i] = zlist[i] - basiswidth[2]*cwlist[i]*.5 + basisheight[2]*chlist[i]
        
        xlistnew[3,i] = xlist[i] - basiswidth[0]*cwlist[i]*.5
        ylistnew[3,i] = ylist[i] - basiswidth[1]*cwlist[i]*.5
        zlistnew[3,i] = zlist[i] - basiswidth[2]*cwlist[i]*.5

    #fig.show()

    #Want to extend the list by a bit on both ends so that the cut completes the chanel
    xlistnew = np.concatenate((np.array([[
        xlist[0]-chanelvectorinit[0]*cwlist[0]*5 + basiswidthinit[0]*cwlist[0]*.5,
        xlist[0]-chanelvectorinit[0]*cwlist[0]*5 + basiswidthinit[0]*cwlist[0]*.5 + basisheightinit[0]*chlist[0],
        xlist[0]-chanelvectorinit[0]*cwlist[0]*5 - basiswidthinit[0]*cwlist[0]*.5 + basisheightinit[0]*chlist[0],
        xlist[0]-chanelvectorinit[0]*cwlist[0]*5 - basiswidthinit[0]*cwlist[0]*.5
    ]]).T,
    xlistnew),axis=1)
    ylistnew = np.concatenate((np.array([[
        ylist[0]-chanelvectorinit[1]*cwlist[0]*5 + basiswidthinit[1]*cwlist[0]*.5,
        ylist[0]-chanelvectorinit[1]*cwlist[0]*5 + basiswidthinit[1]*cwlist[0]*.5 + basisheightinit[1]*chlist[0],
        ylist[0]-chanelvectorinit[1]*cwlist[0]*5 - basiswidthinit[1]*cwlist[0]*.5 + basisheightinit[1]*chlist[0],
        ylist[0]-chanelvectorinit[1]*cwlist[0]*5 - basiswidthinit[1]*cwlist[0]*.5
    ]]).T,
    ylistnew),axis=1)
    zlistnew = np.concatenate((np.array([[
        zlist[0]-chanelvectorinit[2]*cwlist[0]*5 + basiswidthinit[2]*cwlist[0]*.5,
        zlist[0]-chanelvectorinit[2]*cwlist[0]*5 + basiswidthinit[2]*cwlist[0]*.5 + basisheightinit[2]*chlist[0],
        zlist[0]-chanelvectorinit[2]*cwlist[0]*5 - basiswidthinit[2]*cwlist[0]*.5 + basisheightinit[2]*chlist[0],
        zlist[0]-chanelvectorinit[2]*cwlist[0]*5 - basiswidthinit[2]*cwlist[0]*.5
    ]]).T,
    zlistnew),axis=1)

    xlistnew = np.concatenate((
    xlistnew,
    np.array([[xlist[i]+chanelvector[0]*cwlist[i]*.5*5 + basiswidth[0]*cwlist[i]*.5,
        xlist[i]+chanelvector[0]*cwlist[i]*.5*5 + basiswidth[0]*cwlist[i]*.5 + basisheight[0]*chlist[i],
        xlist[i]+chanelvector[0]*cwlist[i]*.5*5 - basiswidth[0]*cwlist[i]*.5 + basisheight[0]*chlist[i],
        xlist[i]+chanelvector[0]*cwlist[i]*.5*5 - basiswidth[0]*cwlist[i]*.5
    ]]).T
    ),axis=1)
    ylistnew = np.concatenate((ylistnew,
    np.array([[
        ylist[i]+chanelvector[1]*cwlist[i]*.5*5 + basiswidth[1]*cwlist[i]*.5,
        ylist[i]+chanelvector[1]*cwlist[i]*.5*5 + basiswidth[1]*cwlist[i]*.5 + basisheight[1]*chlist[i],
        ylist[i]+chanelvector[1]*cwlist[i]*.5*5 - basiswidth[1]*cwlist[i]*.5 + basisheight[1]*chlist[i],
        ylist[i]+chanelvector[1]*cwlist[i]*.5*5 - basiswidth[1]*cwlist[i]*.5
    ]]).T
    ),axis=1)
    zlistnew = np.concatenate((
    zlistnew,
    np.array([[
        zlist[i]+chanelvector[2]*cwlist[i]*.5*5 + basiswidth[2]*cwlist[i]*.5,
        zlist[i]+chanelvector[2]*cwlist[i]*.5*5 + basiswidth[2]*cwlist[i]*.5 + basisheight[2]*chlist[i],
        zlist[i]+chanelvector[2]*cwlist[i]*.5*5 - basiswidth[2]*cwlist[i]*.5 + basisheight[2]*chlist[i],
        zlist[i]+chanelvector[2]*cwlist[i]*.5*5- basiswidth[2]*cwlist[i]*.5
    ]]).T
    ),axis=1)

    return xlistnew, ylistnew,zlistnew

"""#The bean concept is currently as follows:
The temps are calculated assuming box chanels
The bean that best emulates the box is one that is both the same cross sectional area
and the same hydrualic diamemter and the same width. We just have to find the bean
that satisfies all these! We simply solve for all of them being equal with a unique
solution given by adjusting three variables (three unkowns): width, height, and top ellipse vertical axis
The actual bean is normalized to having a width of 1, with the bottom half being the reniform bean
and the top half being an ellipse"""
def ChanelBean(xlistoriginal,rlist,twlist,helicitylist, chlist, cwlist):
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

        

    ylist = np.zeros(len(xlistoriginal))
    ylist[0] = xlistoriginal[0]
    xlist=np.zeros(len(ylist))
    zlist=np.zeros(len(ylist))
    xlist[0]=rlist[0]+twlist[0]
    theta = 0
    for i in range(1,len(ylist)):
        ylist[i] = xlistoriginal[i]
        r=rlist[i]+twlist[i]
        #c=math.tan(helicitylist[i])*2*math.pi*r
        dy=(ylist[i]-ylist[i-1])
        if helicitylist[i]==math.pi/2:
            dcircum = 0
        else:
            dcircum = dy/math.tan(helicitylist[i])
        dtheta = dcircum/r
        theta=theta+dtheta
        xlist[i]=r*math.cos(theta)
        zlist[i]=r*math.sin(theta)
    # Now we hzve the geometry for the sweep curve, or the curve along rlist+twlist
    #We want to make four lists now, one for each corner of the box
    numbeanpoints = 21
    ylistnew = np.zeros((2*numbeanpoints,len(xlistoriginal)))
    xlistnew=np.zeros((2*numbeanpoints,len(ylist)))
    zlistnew =np.zeros((2*numbeanpoints,len(ylist)))
    #fig = plt.figure() # for showing vectors for debugging
    #ax=fig.add_subplot(projection='3d')
    hydraulicdiamlistnew =np.zeros(len(ylist))
    chlistnew =np.zeros(len(ylist))
    for i in range(0,len(ylist)):
        angle = abs(math.atan(xlist[i]/zlist[i]))
        #xlist[i]=xlist[i]+math.sin(angle)*np.sign(xlist[i])*chlist[i]
        #zlist[i]=zlist[i]+math.cos(angle)*np.sign(zlist[i])*chlist[i]
        #Find the angle of rlist (converging angle, divering angle, etc)
        if i%5 == 0: #this is just for debugging
            print(i)
        if i == 0:
            drdy = (rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i])
        elif i == len(xlistoriginal)-1:
            drdy = (rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])
        else:
            drdy = ((rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])+(rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i]))/2
        axialangle = -math.atan(drdy)
        # Now find the basis vectors for the plane with the curve as its normal
        #basis1 = [math.sin(angle)*np.sign(xlist[i]), 0 , math.cos(angle)*np.sign(zlist[i])]
        #basis1*basis2 = 0
        #basis2*[0 1 0] = cos(90-rlistangle)*|basis2|
        #basis2 = [A*math.cos(angle)*np.sign(zlist[i]), cos(math.pi/2-rlistangle)*|basis2|, -A*math.sin(angle)*np.sign(xlist[i])]
        #Set |basis2| =1
        #sqrt(A^2+cos(math.pi/2-rlistangle)^2)=1
        #A^2=1-cos(math.pi/2-rlistangle)^2
        Aaxial=math.sqrt(1-(math.cos(math.pi/2-axialangle)**2))
        #Ahelicity=math.sqrt(1-(math.cos(helicitylist[i])**2))
        basisheight = [Aaxial*math.sin(angle)*np.sign(xlist[i]), math.cos(math.pi/2-axialangle), Aaxial*math.cos(angle)*np.sign(zlist[i])]
        #now with basis height, we want to find a perpendicular vector that is also angled with the helicity
        # for a helicity of 90, the axial (y) component should be zero
        #chaenlcurve vector defined in cylydnircal coordinates : (r,theta,y)
        #-theta/sin(hel)=y/cos(hel) 
        #r/sin(axial)=y/cos(axial)
        #r^2+theta^2+y^2=1
        #ytan(hel)^2+ytan(axial)^2+y^2=1
        #y=math.sqrt(1/(1+math.tan(math.pi/2-helicitylist[i])**2+math.tan(axialangle)**2))
        # x = r*sin(angle) +theta*cos(angle)
        # z = r*cos(angle) - theta*sin(agnle)
        # chanelvector = [math.tan(axialangle)*y*math.sin(angle) - math.tan(helicitylist[i])*y*math.cos(angle)
        # , y, 
        # math.tan(axialangle)*y*math.cos(angle) + math.tan(helicitylist[i])*y*math.sin(angle)]
        #cylyndrilca basis height(r,theta,z) = (cos(axialangle),0,sin(axialangle))
        #cylidrilca chanel vector = [-tan(axialangle)y, -tan(hel)y, y]
        angle = np.arctan2(xlist[i],zlist[i]) # this is so ham but I did it janky before, now I fix
        if helicitylist[i] == math.pi/2:
            y=math.sqrt(1/(1+math.tan(axialangle)**2))
            chanelvector = [-(math.tan(axialangle)*y*math.sin(angle)),
            y, 
            -(math.tan(axialangle)*y*math.cos(angle))]
        else:
            y=math.sqrt(1/(1+(1/math.tan(helicitylist[i]))**2+math.tan(axialangle)**2))
            chanelvector = [-(math.tan(axialangle)*y*math.sin(angle) +(1/math.tan(helicitylist[i]))*y*math.cos(angle)),
            y, 
            -(math.tan(axialangle)*y*math.cos(angle) - (1/math.tan(helicitylist[i]))*y*math.sin(angle))]
        
        

        # the perpendicular vector is just chanelvector cross basisheight
        basiswidth = np.cross(chanelvector,basisheight)
        if i == 0: # this is for the exhaust-side extension (to ensure complete cutting)
            chanelvectorinit = chanelvector
            basisheightinit = basisheight
            basiswidthinit = basiswidth
        #Helps verify vectors point where they should
        #ax.plot([xlist[i],xlist[i] + basiswidth[0]*.01],[ylist[i],ylist[i] + basiswidth[1]*.01],
        #                [zlist[i],zlist[i] + basiswidth[2]*.01],'r')
        #ax.plot([xlist[i],xlist[i] + chanelvector[0]*.01],[ylist[i],ylist[i] + chanelvector[1]*.01],
        #                [zlist[i],zlist[i] + chanelvector[2]*.01],'g')
        #ax.plot([xlist[i],xlist[i] + basisheight[0]*.01],[ylist[i],ylist[i] + basisheight[1]*.01],
        #                [zlist[i],zlist[i] + basisheight[2]*.01],'b')

        
        #basiswidth = [Ahelicity*math.cos(angle)*np.sign(zlist[i]), math.cos(helicitylist[i]), -Ahelicity*math.sin(angle)*np.sign(xlist[i])]
        
        #numbeanpoints = xlistnew.shape[0]/4 # a bean points defines 2 actual points, one in each half (pos and neg)
        beanpoints = np.linspace(-.45,.45,numbeanpoints)
        beanpoints = np.vstack((beanpoints,np.zeros((2,numbeanpoints))))
        beanpoints[1,:] = ellipsefunc(beanpoints[0,:],axisy = .35,axisx = .5, sign=1)+.5
        beanpoints[2,:] =  ellipsefunc(beanpoints[0,:],axisy = .4,axisx = .5, sign=-1)+.5 + .175*(np.cos(beanpoints[0,:]*math.pi*2)+1)
        beanreduce = np.min(beanpoints[2,:])
        beanpoints[1,:] = beanpoints[1,:] - beanreduce
        beanpoints[2,:]= beanpoints[2,:] - beanreduce
        beanmultiply = np.max(beanpoints[1,:])
        beanpoints[1,:] = beanpoints[1,:]/beanmultiply
        beanpoints[2,:]= beanpoints[2,:] /beanmultiply
        beanpoints[1,:] = beanpoints[1,:]*chlist[i]
        beanpoints[2,:] = beanpoints[2,:]*chlist[i]
        beanpoints[0,:] = beanpoints[0,:]*cwlist[i]
        areadesired = chlist[i]*cwlist[i]
        perimdesired = 2*chlist[i]+2*cwlist[i]
        axisyval = .2
        chguess = chlist[i]
        areacurrent = area(beanpoints)
        perim = perimeter(beanpoints)
        tol=.0001
        while abs(areacurrent-areadesired)/areadesired>tol:
            beanpoints = np.linspace(-.45,.45,numbeanpoints)
            beanpoints = np.vstack((beanpoints,np.zeros((2,numbeanpoints))))
            chguess = chguess- (areacurrent-areadesired)/areadesired*chguess
            beanpoints[1,:] = ellipsefunc(beanpoints[0,:],axisy = .35,axisx = .5, sign=1)+.5
            beanpoints[2,:] =  ellipsefunc(beanpoints[0,:],axisy = .4,axisx = .5, sign=-1)+.5 + .175*(np.cos(beanpoints[0,:]*math.pi*2)+1)
            beanreduce = np.min(beanpoints[2,:])
            beanpoints[1,:] = beanpoints[1,:] - beanreduce
            beanpoints[2,:]= beanpoints[2,:] - beanreduce
            beanmultiply = np.max(beanpoints[1,:])
            beanpoints[1,:] = beanpoints[1,:]/beanmultiply
            beanpoints[2,:]= beanpoints[2,:] /beanmultiply
            beanpoints = np.vstack((beanpoints,np.zeros((2,numbeanpoints))))
            beanpoints[1,:] = beanpoints[1,:]*chguess
            beanpoints[2,:] = beanpoints[2,:]*chguess
            beanpoints[0,:] = beanpoints[0,:]*cwlist[i]
            chguess = chguess- (areacurrent-areadesired)/areadesired*chguess
            areacurrent = area(beanpoints)
            perim = perimeter(beanpoints)
        #plt.plot(beanpoints[0,:],beanpoints[1,:])
        #plt.plot(beanpoints[0,:],beanpoints[2,:])
        #beanpoints[1,:] = -np.min(beanpoints[2,:])+beanpoints[1,:] #currently bottom is at -.5*ch, gotta shift it up
        #beanpoints[2,:] = -np.min(beanpoints[2,:]) + beanpoints[2,:]
        #plt.plot(beanpoints[0,:],beanpoints[1,:])
        #plt.plot(beanpoints[0,:],beanpoints[2,:])
        #plt.show()
        hydraulicdiamlistnew[i] = 4*areacurrent/perim
        chlistnew[i] = np.max(beanpoints[1,:])
        if i == 0: # saving for future extension
            beanpointsinit = beanpoints
        for beanindex in range(0,numbeanpoints):
            
            xlistnew[beanindex,i] = xlist[i] + basisheight[0]*beanpoints[1,beanindex] + basiswidth[0]*beanpoints[0,beanindex]
            ylistnew[beanindex,i] = ylist[i] + basisheight[1]*beanpoints[1,beanindex] + basiswidth[1]*beanpoints[0,beanindex]
            zlistnew[beanindex,i] = zlist[i] + basisheight[2]*beanpoints[1,beanindex] + basiswidth[2]*beanpoints[0,beanindex]

            xlistnew[numbeanpoints*2-beanindex-1,i] = xlist[i] + basisheight[0]*beanpoints[2,beanindex] + basiswidth[0]*beanpoints[0,beanindex]
            ylistnew[numbeanpoints*2-beanindex-1,i] = ylist[i] + basisheight[1]*beanpoints[2,beanindex] + basiswidth[1]*beanpoints[0,beanindex]
            zlistnew[numbeanpoints*2-beanindex-1,i] = zlist[i] + basisheight[2]*beanpoints[2,beanindex] + basiswidth[2]*beanpoints[0,beanindex]
        print(i)



    #fig.show()

    #Want to extend the list by a bit on both ends so that the cut completes the chanel
    precatarray = np.zeros((numbeanpoints*2,1))
    for beanindex in range(0,numbeanpoints):
            
            precatarray[beanindex,0] = xlist[0] -chanelvectorinit[0]*cwlist[0]+ basisheightinit[0]*beanpointsinit[1,beanindex] + basiswidthinit[0]*beanpointsinit[0,beanindex]
            precatarray[numbeanpoints*2-beanindex-1,0] = xlist[0] -chanelvectorinit[0]*cwlist[0]+ basisheightinit[0]*beanpointsinit[2,beanindex] + basiswidthinit[0]*beanpointsinit[0,beanindex]


    xlistnew = np.concatenate((precatarray,
    xlistnew),axis=1)

    for beanindex in range(0,numbeanpoints):
            
            precatarray[beanindex,0] = ylist[0] -chanelvectorinit[1]*cwlist[0]+ basisheightinit[1]*beanpointsinit[1,beanindex] + basiswidthinit[1]*beanpointsinit[0,beanindex]
            precatarray[numbeanpoints*2-beanindex-1,0] = ylist[0] -chanelvectorinit[1]*cwlist[0]+ basisheightinit[1]*beanpointsinit[2,beanindex] + basiswidthinit[1]*beanpointsinit[0,beanindex]

    ylistnew = np.concatenate((precatarray,
    ylistnew),axis=1)
    
    for beanindex in range(0,numbeanpoints):
            
            precatarray[beanindex,0] = zlist[0] -chanelvectorinit[2]*cwlist[0]+ basisheightinit[2]*beanpointsinit[1,beanindex] + basiswidthinit[2]*beanpointsinit[0,beanindex]
            precatarray[numbeanpoints*2-beanindex-1,0] = zlist[0] -chanelvectorinit[2]*cwlist[0]+ basisheightinit[2]*beanpointsinit[2,beanindex] + basiswidthinit[2]*beanpointsinit[0,beanindex]

    zlistnew = np.concatenate((precatarray,
    zlistnew),axis=1)

    postcatarray = np.zeros((numbeanpoints*2,1))
    for beanindex in range(0,numbeanpoints):
            
            postcatarray[beanindex,0] = xlist[i] +chanelvector[0]*cwlist[i]+ basisheight[0]*beanpoints[1,beanindex] + basiswidth[0]*beanpoints[0,beanindex]
            postcatarray[numbeanpoints*2-beanindex-1,0] = xlist[i] +chanelvector[0]*cwlist[i]+ basisheight[0]*beanpoints[2,beanindex] + basiswidth[0]*beanpoints[0,beanindex]

    xlistnew = np.concatenate((
    xlistnew,
    postcatarray
    ),axis=1)
    
    for beanindex in range(0,numbeanpoints):
            
            postcatarray[beanindex,0] = ylist[i] +chanelvector[1]*cwlist[i]+ basisheight[1]*beanpoints[1,beanindex] + basiswidth[1]*beanpoints[0,beanindex]
            postcatarray[numbeanpoints*2-beanindex-1,0] = ylist[i] +chanelvector[1]*cwlist[i]+ basisheight[1]*beanpoints[2,beanindex] + basiswidth[1]*beanpoints[0,beanindex]

    ylistnew = np.concatenate((
    ylistnew,
    postcatarray
    ),axis=1)

    for beanindex in range(0,numbeanpoints):
            
            postcatarray[beanindex,0] = zlist[i] +chanelvector[2]*cwlist[i]+ basisheight[2]*beanpoints[1,beanindex] + basiswidth[2]*beanpoints[0,beanindex]
            postcatarray[numbeanpoints*2-beanindex-1,0] = zlist[i] +chanelvector[2]*cwlist[i]+ basisheight[2]*beanpoints[2,beanindex] + basiswidth[2]*beanpoints[0,beanindex]

    zlistnew = np.concatenate((
    zlistnew,
    postcatarray
    ),axis=1)


    return xlistnew, ylistnew,zlistnew, hydraulicdiamlistnew, chlistnew

def rlistExtender(xlist,rlist,ewlist):
    externalrlist = np.zeros(len(rlist))
    newxlist = np.zeros(len(xlist))
    normallist = np.arctan2(np.diff(rlist),np.diff(xlist))+math.pi/2
    normallist = np.concatenate(([normallist[0]],normallist))
    externalrlist = rlist+np.sin(normallist)*ewlist
    newxlist = xlist+np.cos(normallist)*ewlist
    return newxlist,externalrlist


    



def ChanelGuidingCurve_Height(xlistoriginal,rlist,twlist,helicitylist, chlist):
    ylist = np.zeros(len(xlistoriginal))
    xlist=np.zeros(len(ylist))
    zlist=np.zeros(len(ylist))
    xlist[0]=rlist[0]+twlist[0]
    theta = 0
    for i in range(1,len(ylist)):
        ylist[i] = xlistoriginal[i]
        r=rlist[i]+twlist[i]
        c=math.tan(helicitylist[i])*2*math.pi*r
        dy=(ylist[i]-ylist[i-1])
        dcircum = dy/math.tan(helicitylist[i])
        dtheta = dcircum/r
        theta=theta+dtheta
        xlist[i]=r*math.cos(theta)
        zlist[i]=r*math.sin(theta)
    # Now we hzve the geometry for the sweep curve, or the curve along rlist+twlist
    #We want to extend the curve along r by ch
    for i in range(0,len(ylist)):
        angle = abs(math.atan(xlist[i]/zlist[i]))
        #xlist[i]=xlist[i]+math.sin(angle)*np.sign(xlist[i])*chlist[i]
        #zlist[i]=zlist[i]+math.cos(angle)*np.sign(zlist[i])*chlist[i]
        #Find the angle of rlist (converging angle, divering angle, etc)
        if i%50 == 0:
            i=i
        if i == 0:
            drdy = (rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i])
        elif i == len(xlistoriginal)-1:
            drdy = (rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])
        else:
            drdy = ((rlist[i]-rlist[i-1])/(xlistoriginal[i]-xlistoriginal[i-1])+(rlist[i+1]-rlist[i])/(xlistoriginal[i+1]-xlistoriginal[i]))/2
        rlistangle = abs(math.atan(drdy))
        # Now find the basis vectors for the plane with the curve as its normal
        #basis1 = [math.sin(angle)*np.sign(xlist[i]), 0 , math.cos(angle)*np.sign(zlist[i])]
        #basis1*basis2 = 0
        #basis2*[0 1 0] = cos(90-rlistangle)*|basis2|
        #basis2 = [A*math.cos(angle)*np.sign(zlist[i]), cos(math.pi/2-rlistangle)*|basis2|, -A*math.sin(angle)*np.sign(xlist[i])]
        #Set |basis2| =1
        #sqrt(A^2+cos(math.pi/2-rlistangle)^2)=1
        #A^2=1-cos(math.pi/2-rlistangle)^2
        A=math.sqrt(1-(math.cos(math.pi/2-rlistangle)**2))
        basis2 = [A*math.sin(angle)*np.sign(xlist[i]), math.cos(math.pi/2-rlistangle), A*math.cos(angle)*np.sign(zlist[i])]
        xlist[i] = xlist[i] + basis2[0]*chlist[i]
        ylist[i] = ylist[i] + basis2[1]*chlist[i]
        zlist[i] = zlist[i] + basis2[2]*chlist[i]
    return xlist, ylist, zlist

def ChanelGuidingCurve_Width(xlistoriginal,rlist,twlist,helicitylist, cwlist):
    ylist = np.zeros(len(xlistoriginal))
    xlist=np.zeros(len(ylist))
    zlist=np.zeros(len(ylist))
    xlist[0]=rlist[0]+twlist[0]
    theta = 0
    for i in range(1,len(ylist)):
        ylist[i]=xlistoriginal[i]
        r=rlist[i]+twlist[i]
        c=math.tan(helicitylist[i])*2*math.pi*r
        dy=(ylist[i]-ylist[i-1])
        dcircum = dy/math.tan(helicitylist[i])
        dtheta = dcircum/r
        theta=theta+dtheta
        xlist[i]=r*math.cos(theta)
        zlist[i]=r*math.sin(theta)
    for i in range(0,len(ylist)):
        # start with x and z in normal plane of sweep curve
        angle = abs(math.atan(xlist[i]/zlist[i]))
        #xnormal=xlist[i]+math.cos(angle)*np.sign(zlist[i])*cwlist[i]
        #znormal=zlist[i]-math.sin(angle)*np.sign(xlist[i])*cwlist[i]
        # Now find the basis vectors for the plane with the curve as its normal
        #basis1 = [math.sin(angle)*np.sign(xlist[i]), 0 , math.cos(angle)*np.sign(zlist[i])]
        #basis1*basis2 = 0
        #basis2*[0 1 0] = cos(90-helicity)*|basis2|
        #basis2 = [A*math.cos(angle)*np.sign(zlist[i]), cos(math.pi/2-helicitylist[i])*|basis2|, -A*math.sin(angle)*np.sign(xlist[i])]
        #Set |basis2| =1
        #sqrt(A^2+cos(math.pi/2-helicitylist[i])^2)=1
        #A^2=1-cos(math.pi/2-helicitylist[i])^2
        A=math.sqrt(1-(math.cos(math.pi/2-helicitylist[i])**2))
        basis2 = [A*math.cos(angle)*np.sign(zlist[i]), math.cos(math.pi/2-helicitylist[i]), -A*math.sin(angle)*np.sign(xlist[i])]
        xlist[i] = xlist[i] + basis2[0]*cwlist[i]
        ylist[i] = ylist[i] + basis2[1]*cwlist[i]
        zlist[i] = zlist[i] + basis2[2]*cwlist[i]
    return xlist, ylist, zlist

def ChanelSweepCurve(ylist,rlist,twlist,helicitylist): #outdated
    xlist=np.zeros(len(ylist))
    zlist=np.zeros(len(ylist))
    xlist[0]=rlist[0]+twlist[0]
    theta = 0
    for i in range(1,len(ylist)):
        r=rlist[i]+twlist[i]
        c=math.tan(helicitylist[i])*2*math.pi*r
        dy=(ylist[i]-ylist[i-1])
        dcircum = dy/math.tan(helicitylist[i])
        dtheta = dcircum/r
        theta=theta+dtheta
        xlist[i]=r*math.cos(theta)
        zlist[i]=r*math.sin(theta)
    return xlist, ylist, zlist


