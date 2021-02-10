# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:19:15 2016

@author: coulais
"""
import matplotlib as plt
import cv2
import csv
from scipy import ndimage,sparse
import os
import matplotlib as mpl
import numpy as np
import sys
from scipy.optimize import curve_fit
from scipy.cluster.vq import kmeans, vq
import networkx as nx
from matplotlib.mlab import find


def measure1(imgcrop,threshold=200,threshold2=255,ax1=None,ax2=None):
    (tmp,ThresholdedImage)= cv2.threshold(imgcrop,threshold,threshold2,cv2.THRESH_BINARY_INV)

    tmp,contours, hierarchy = cv2.findContours(255-ThresholdedImage,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    del tmp
    s=np.array([len(c) for c in contours])
    #print(s)
    x,y,w,h = cv2.boundingRect(contours[s.argmax()])
    BBox=cv2.boundingRect(contours[s.argmax()])

    Area=(255-ThresholdedImage[y:y+h,x:x+w]).sum()/(255)
    if ax2 is not None:
        ax2.imshow(255-ThresholdedImage[y:y+h,x:x+w])
    if ax1 is not None:
            ax1.plot([x,x+w,x+w,x,x],[y,y,y+h,y+h,y],color="w")
    return Area, w,BBox

def measure2(img,BBox,dx=100,dy=2,ax1=None,ax2=None):
    x,y,w,h=BBox

    imgcrop2=img[y+dy:y+h-dy,x:x+w+dx]
    profile=imgcrop2.sum(axis=0)/(255.*np.size(imgcrop2[:,0]))
    x=np.array(range(len(profile)))
    #a=0.3;b=0.22;c=1;d=120.
    #print profile
    popt, pcov = curve_fit(step_atan, x,profile,maxfev = 10000)
    a,b,c,d=popt
    #print popt
    if ax2 is not None:
        ax2.plot(profile,color="b")
        ax2.plot(step_atan(x,a,b,c,d),color="r")
    return d/c
    '''
    Area=(255-ThresholdedImage[y:y+h,x:x+w]).sum()/(255)
    if ax2 is not None:
        ax2.imshow(255-ThresholdedImage[y:y+h,x:x+w])
    if ax1 is not None:
            ax1.plot([x,x+w,x+w,x,x],[y,y,y+h,y+h,y],color="w")
    return Area, w
    '''
def measure3(img,BBox,dx=100,dy=2,ax1=None,ax2=None):
    x,y,w,h=BBox

    imgcrop2=img#[y+dy:y+h-dy,x:x+w+dx]
    profile=imgcrop2.sum(axis=0)/(255.*np.size(imgcrop2[:,0]))
    x=np.array(range(len(profile)))
    #a=0.3;b=0.22;c=1;d=120.
    #print profile
    #plt.plot(profile)
    popt, pcov = curve_fit(gaussian, x,profile,maxfev = 10000,p0=[0.1,0.7,1e-2,100])#,method='trf',bounds=([0.,0.5,1e-3,0],[0.,1,1e-1,200]))
    a,b,c,d=popt
    #print popt
    if ax2 is not None:
        ax2.plot(profile,color="b")
        ax2.plot(gaussian(x,a,b,c,d),color="r")
    return d

def measure4(img,BBox,dx=100,dy=2,ths=0.35,ax1=None,ax2=None):
    x,y,w,h=BBox

    imgcrop2=img#[y+dy:y+h-dy,x:x+w+dx]
    profile=imgcrop2.sum(axis=0)/(255.*np.size(imgcrop2[:,0]))
    x=np.array(range(len(profile)))
    #find clusters and keep only the bigger
    idx=np.array(np.where(profile>ths)[0])
    didx=np.diff(idx)
    clus=[i for i in range(len(didx)) if didx[i]>1]
    #print(idx,didx,clus)
    clusters=[]
    if len(clus)== 0:
        clusters+=[idx]
        im=0
    else:
        clusters+=[idx[:clus[0]+1]]
        for cc in range(1,len(clus)-1):
            clusters+=[idx[clus[cc]+1:clus[cc+1]+1]]
        clusters+=[idx[clus[-1]+1:]]
        im=np.array([len(c) for c in clusters]).argmax()
        #print(clusters[im])
    profile2=np.zeros(len(x))
    profile2[clusters[im]]=profile[clusters[im]]

    # try:
    #     popt, pcov = curve_fit(gaussian, x,profile2,p0=[0.1,0.8,2e-3,300],method='trf',bounds=([0.02,0.5,1e-3,0],[0.2,1,1e-2,700]))
    #     a,b,c,d=popt
    #except:
    #     d=np.NaN
    d=np.NaN
    d2=np.sum(profile2*x)/profile2.sum()
    #print popt
    if ax2 is not None:
        ax2.plot(profile,color="b")
        ax2.plot(profile2,color="LightBlue")
        ax2.plot(gaussian(x,a,b,c,d),color="r")
        ax2.vlines(x=d2,ymin=0,ymax=1,color="r")
    return d,d2

def step_atan(x,a,b,c,d):
    #a=0.3;b=0.22;c=1;d=120.
    #a,b,c,d=params
    return a-b*np.arctan(x*c-d)

def gaussian(x,a,b,c,d):
    #a=0.3;b=0.22;c=1;d=120.
    #a,b,c,d=params
    return a+b*np.exp(-c*(x-d)**2)

def fitellipses(ThresholdedImage,lim,BBox=[0,0,0,0]):
    x0,deltax,y0,deltay=BBox

    ######## Find contours ####
    im2, contours, hierarchy = cv2.findContours(ThresholdedImage,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    ellipses = []
    #fitEllipse needs more than 4 datapoints to fit an ellipse
    for ii in np.arange(0,len(contours)):
        if len(contours[ii])>lim:#Fit an ellipse for all contours with length>4
            ellips = cv2.fitEllipse(contours[ii])
            ellipse =  ((ellips[0][0]+x0,ellips[0][1]+y0),(ellips[1]),ellips[2])
            ellipses.append(ellipse)#ellipse = ((x0,y0),(a,b),angle)
    return (ellipses)

def save_ellipses(jj,ellipses,img,BBox,minArea,maxArea,maxaspect,ellipses_out=[],bright=255):
    x0,deltax,y0,deltay=BBox
    num = 0
    for kk in np.arange(0,len(ellipses)):
        #area of ellipse = pi* a *b
        area_ellipse = ellipses[kk][1][1]*ellipses[kk][1][0]
        major_radius = ellipses[kk][1][1]
        minor_radius = ellipses[kk][1][0]
        aspectratio=major_radius/minor_radius
        if (area_ellipse < maxArea and area_ellipse > minArea and aspectratio < maxaspect):
            num =num+1
            ellipse = ellipses[kk]
            #print ellipse
            ellipse2 = ((np.copy(ellipse[0][0]-x0),np.copy(ellipse[0][1]-y0)),(np.copy(ellipse[1][0]+10),np.copy(ellipse[1][1]+10)),np.copy(ellipse[2]))
            #cv2.ellipse(CLONE,ellipse,(bright,bright,bright),-1)#draw this ellips to check
            cv2.ellipse(img,ellipse2,(bright,bright,bright),-1)#draw this ellips to check
            ellipses_out+=[ellipse]
            #text_file.writelines([str(jj)+'\t',str(num)+'\t',str(ellipses[kk][0][0])+'\t',str(ellipses[kk][0][1])+'\t',str(ellipses[kk][1][0])+'\t',str(ellipses[kk][1][1])+'\t',str(ellipses[kk][2])+'\n'])
    return (ellipses_out)

    #cv2.imwrite(pathchecks+file_out,CLONE)

def plot_ellipses(ax,ellipses,CLONE,pathchecks,file_out):
    for kk in np.arange(0,len(ellipses)):
        #area of ellipse = pi* a *b
        ellipse = ellipses[kk]
        cv2.ellipse(ax,CLONE,ellipse,(0,255,0),2)#draw this ellips to check

    cv2.imwrite(pathchecks+file_out,CLONE)

def tracking(tpos,mem=1,delta=1,display=False,usesparse=True):
    #Tracking program:
    # you need to input a array tpos =[t X Y], with one particle for each line
    # the ouput is a list containining distinct [t X Y] arrays for each tracked particle

    # options include:
    #     the maximal displacement between two frame delta
    #     the memory of the tracking
    #     display step number and number of particles
    #     the use of sparse matrix (a bit slower but allows to track more particles and more time steps)

    ######## typical use #########
    #data=np.loadtxt('/home/science/geek/4ellipticity.txt',skiprows=1)
    #X=data[:,2];Y=data[:,3];t=data[:,0]
    #tpos=data[:,[0,2,3]]

    #P=tracking(tpos)
    #########################

    #Prerequisites:
    ### YOU NEED TO INSTALL THE Python-networkx library ###

    ## Corentin Coulais, 23/04/2014, last updated 19/06/2014

   #mem=1;delta=100
   # --- Definitions
    
    t = tpos[:,0]
    pos = tpos[:,2:4:];
    extra = tpos[:,4::];
    n, dim = pos.shape
    print(n, dim)

    # --- Preparation

    t_ = np.unique(t);

    if usesparse == False:
        A=np.zeros((n,n))
    else:
        A = sparse.csr_matrix((n,n))
        #A = sparse.lil_matrix((n,n))

    for i in range(len(t_)-mem):

        #Get the particules at t and in ]t;t+mem]
        I = find(t==t_[i]);
        J = find((t>t_[i])*(t<=t_[i]+mem));

        # Exclude already linked
        if usesparse == False:
            J = J[~A[:,J].any(0)]
        else:#not so elegant hack for sparse matrices, a bit slow I think
            JJ=[]
            for j in J:
                if len(A.getcol(j).nonzero()[0])==0:
                    JJ+=[j]
            J=np.array(JJ)


        ### Create distance and time steps matrix ###
        [J_, I_] = np.meshgrid(J,I);
        K_ = (t[J_] - t[I_]).reshape((len(I),len(J)))
        D=((pos[I_,:] - pos[J_,:])**2).sum(2); D=np.sqrt(D)
        D=D.reshape((len(I),len(J)))
        ### compute distance### if you divide out by K_ thn you put a treshold on velocities
        V = D#/K_;
        ### Apply displacement threshold ###
        #L=np.array(V<delta,int);
        L=np.array(np.floor((np.array(V<delta,int)+np.array(K_!=0,int))/2),int)
        #### Keep the closests in time ###
        prod=L*K_
        L=prod==np.ones((len(I),1))*prod.min(0)    #L = (L*K_==np.ones(len(I),1)*min(L*K_,[],1))*1;

        ### Keep the closest in distance ###
        prod=L*D
        #print prod
        M=prod.argmin(0)
        M[prod.min(0)>=delta]=-1

        ### Build adjacency matrix (associate the index of linked particles) ###
        for j in range(len(M)):
            if M[j] != -1:
                A[I[M[j]], J[j]] = 1;

        ### optional display ###
        if display:
            print('step %d, %d tracked particles' % (i,len(M)))

    #### symetrise A ####
    A=A+A.transpose()

    ## attempt to use sparse matrix (NOT USED YET)##
    #B=scipy.sparse.csr_matrix(A)

    ### create a graph and use the networkx toolbox ###
    D=nx.DiGraph(A)
    #del A
    ### convert the graph into a list of connected clusters ###
    #C=nx.connected.connected_components(D.to_undirected())
    C=sorted(nx.connected_components(D.to_undirected()))
    #del D
    ### loop over them and build a list of tracked particle ###
    P = []
    for i in range(len(C)):
        #ix=t[C[i]].argsort()
        #idx=np.array(C[i])[ix]
        L=np.array(list(C[i]))
        ix=t[L].argsort()
        idx=L[ix]
        if len(L)==1:### small hack when there is only one particle
            tmp=np.hstack((np.transpose(t[idx]),pos[idx[0],:],extra[idx[0],:]))
        elif len(L) >1:
            tmp=np.transpose(np.vstack((t[idx],np.transpose(pos[idx,:]),np.transpose(extra[idx,:]))))
        P+=[tmp]

#del C
    del C
    return (P)

def ellipse2point(ellipse,theta,a=100):
    dx=a*np.cos(theta)
    dy=a*np.sin(theta)
    th=ellipse[2]#*np.pi/180
    x2=ellipse[0][0]+dx*np.sin(th)+dy*np.cos(th)
    y2=ellipse[0][1]+-dx*np.cos(th)+dy*np.sin(th)
    return x2,y2

def ellipse2square_old(ellipse,squaresize=50,inte=True):
    phi=ellipse[2]
    pt1=ellipse2point(ellipse,0+phi,a=squaresize)
    pt2=ellipse2point(ellipse,np.pi/2+phi,a=squaresize)
    pt3=ellipse2point(ellipse,np.pi+phi,a=squaresize)
    pt4=ellipse2point(ellipse,-np.pi/2+phi,a=squaresize)
    box=((pt1[0],pt1[1]), (pt2[0], pt2[1]), (pt3[0], pt3[1]), (pt4[0], pt4[1]))
    if inte is True: box = np.int0(box)
    return box

def ellipse2square(ellipse,squaresize=50,inte=True):
    phi=ellipse[2]
    pt1=[ellipse[0][0]+squaresize*np.cos(phi),ellipse[0][1]+squaresize*np.sin(phi)]
    pt2=[ellipse[0][0]+squaresize*np.cos(np.pi/2+phi),ellipse[0][1]+squaresize*np.sin(np.pi/2+phi)]
    pt3=[ellipse[0][0]+squaresize*np.cos(np.pi+phi),ellipse[0][1]+squaresize*np.sin(np.pi+phi)]
    pt4=[ellipse[0][0]+squaresize*np.cos(-np.pi/2+phi),ellipse[0][1]+squaresize*np.sin(-np.pi/2+phi)]
    box=((pt1[0],pt1[1]), (pt2[0], pt2[1]), (pt3[0], pt3[1]), (pt4[0], pt4[1]))
    if inte is True: box = np.int0(box)
    return box

