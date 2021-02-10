# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 19:26:12 2014

@author: coulais
"""

from pylab import ginput
import numpy as np
#from matplotlib.mlab import find

def plotellipse(ax,a,b,phi,x0,y0,linestyle='-',linewidth=1,color='k'):
        ######## plot #######
        x=np.linspace(0,2*np.pi,50)    
        X=a*np.cos(x)*np.cos(phi)+b*np.sin(x)*np.sin(phi)+x0
        Y=b*np.sin(x)*np.cos(phi)-a*np.cos(x)*np.sin(phi)+y0
        ax.plot(X,Y,linestyle=linestyle,linewidth=linewidth,color=color)
        
def plotsquare(ax,a,b,phi,x0,y0,linestyle='-',linewidth=1,color='k'):
        ######## plot #######
        x=np.linspace(0,2*np.pi,4)    
        X=a*np.cos(x)*np.cos(phi)+b*np.sin(x)*np.sin(phi)+x0
        Y=b*np.sin(x)*np.cos(phi)-a*np.cos(x)*np.sin(phi)+y0
        ax.plot(X,Y,linestyle=linestyle,linewidth=linewidth,color=color)

def pick(p,X,Y,N=1):
    print("Please click")

    x=ginput(N)
    print("point",x)
    if N==1:
        ix=np.argmin((X-x[0][0])**2+(Y-x[0][1])**2)
        print("selected particle:",p[ix])
    else:
        ix=np.zeros(N)
        for i in range(N):
            ix[i]=np.argmin((X-x[i][0])**2+(Y-x[i][1])**2)
        ix=np.array(ix,int)
    pout=p[ix]
    print("selected particles:",pout)
    return (pout,ix)
    
def pickrectangle(ax,p,X,Y):
    print("Please click")

    x=ginput(2)
    x1=x[0][0];y1=x[0][1]
    x2=x[1][0];y2=x[1][1]
    ax.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1],'g')
    cond1=np.double((X>=min([x1,x2]))==True)
    cond2=np.double((X<=max([x1,x2]))==True)
    cond3=np.double((Y>=min([y1,y2]))==True)
    cond4=np.double((Y<=max([y1,y2]))==True)
    ix=find(cond1+cond2 +cond3 + cond4==4)
    
    ix=np.array(ix,int)
    pout=p[ix]
    print("selected particles:",pout)
    return (pout,ix)