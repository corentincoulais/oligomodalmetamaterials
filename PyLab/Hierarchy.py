# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:15:30 2015

@author: coulais
"""

import time 
import os, traceback
import numpy as np
from scipy.optimize import minimize
import pickle# as pickle
import matplotlib.pyplot as plt
import cv2
import csv
from scipy import ndimage
import os
import matplotlib as mpl
import numpy as np
import sys
import time
from matplotlib.patches import Polygon,Ellipse
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
from matplotlib.patches import PathPatch
#reload(sys)
#sys.setdefaultencoding("UTF8")
import os
from scipy.optimize import curve_fit
import shutil
from scipy.interpolate import interp1d
from sys import platform as _platform
from matplotlib.mlab import find
from math import isnan


#path_sys={"darwin":"/Volumes/data/AMOLF/groups/hecke-group/Corentin/rawdata",
#          "win32":"//storage01//data//AMOLF//groups//hecke-group//Corentin//rawdata",
#          "linux":"/run/user/1000/gvfs/smb-share:server=storage01,share=data/AMOLF/groups/hecke-group///Corentin//rawdata"}
path_sys={"darwin":"/Users/coulais/science/Hierarchical",
          "win32":"//storage01//data//AMOLF//groups//hecke-group//Corentin//rawdata",
          "linux":"/run/user/1000/gvfs/smb-share:server=storage01,share=data/AMOLF/groups/hecke-group///Corentin//rawdata"}


sys.path.append('/Users/coulais/science/PyLab/')
from GraphicTools import fitellipses,save_ellipses,plot_ellipses, tracking,ellipse2square,ellipse2point
#sys.path.append(path_sys[_platform]+'/PyLocal')
from Data import Database,Part,Time,Prop,Link
#from GraphicTools import fitellipses,save_ellipses,plot_ellipses, tracking,ellipse2square,ellipse2point
sys.path.append('/Users/coulais/science/MyPyBox/')
from Graphics import Tools


def jumpcorrector(datarray,tolerance,tolerance2=30):
    #switch=True
    #jump=0
    jump=180
    for j in range(len(datarray)-1):
        diffj=np.NaN
        if isnan(datarray[j]): continue
        incr=1;diffj=datarray[j]-datarray[j+incr]
        if np.abs(diffj) > tolerance:
            jump=np.rint(diffj/180)*180
            #print(j,incr,diffj,jump)
            datarray[j+incr]+=jump#*np.sign(diffj)
    
    #datarray[np.hstack([np.abs(np.diff(datarray))>tolerance2,False])]=np.NaN
    return datarray
    
def jumpcorrectorold2(datarray,tolerance,tolerance2=30):
    #switch=True
    #jump=0
    jump=180
    for j in range(len(datarray)-1):
        diffj=np.NaN
        if isnan(datarray[j]): continue
        incr=1;diffj=np.abs(datarray[j]-datarray[j+incr])
        while isnan(diffj) is True and incr < np.min([20,len(datarray)-j]):
            incr+=1            
            diffj=np.abs(datarray[j]-datarray[j+incr])                    
        if diffj > tolerance:
            
            jump=np.rint(diffj/180)*180
            print(j,incr,diffj,jump)
            #datarray[j+incr]+=jump*np.sign(datarray[j]-datarray[j+incr])
            datarray[j+incr]+=jump*np.sign(diffj)
    
    #datarray[np.hstack([np.abs(np.diff(datarray))>tolerance2,False])]=np.NaN
    return datarray
    
def jumpcorrectorold(datarray,tolerance,tolerance2=30):
    switch=True
    jump=0
    for j in range(len(datarray)-1):
        if np.abs(datarray[j]-datarray[j+1]) > tolerance: 
            if switch==True:
                jump=180#np.abs(datarray[j]-datarray[j+1])
                #switch=False
            datarray[j+1]+=jump*np.sign(datarray[j]-datarray[j+1])
    
    #idx=        
    datarray[np.hstack([np.abs(np.diff(datarray))>tolerance2,False])]=np.NaN
    #datarray[np.abs(np.diff(datarray))>tolerance2]=np.NaN
    return datarray

### subclasses for Focus ###
class Infos(object):
    pass

### subclasses for Focus ###
class Movie(object):
    pass

class Timeseries(object):
    pass

class Imagedata(object):
    pass

class Tmp(object):
    pass

class ROI(object):
    pass

#class Link(object):
#    pass


### big class Focus ###
class Focus(object):
    #path='/home/coulais/science/holographicelasticity/'
    #path=os.getcwd()+'/'
    path=path_sys[_platform]+'/'
    print(path)
    
    #path='/run/user/1000/gvfs/smb-share:server=storage01,share=data/AMOLF/groups/hecke-group//Corentin//rawdata//'
    name=''
    sub=np.NaN
    def __init__(self,name,sub,display=False):
        self.path=path_sys[_platform]+'/'#"/Volumes/data/AMOLF/groups/hecke-group/Corentin/rawdata/"        
        self.name=name
        self.sub=sub      
        self.DB=Database.DB()
        self.DB.connect(self.name)
        
        self.infos=Infos()
        self.infos.file=self.path + self.name + '//Files//infos' + str(self.sub) + '.p'
        self.load(dat='infos')
        self.infos.file=self.path + self.name + '//Files//infos' + str(self.sub) + '.p'
        #self.infos.ROI=ROI()
        
        self.timeseries=Timeseries()
        self.timeseries.file=self.path + self.name + '//Files//timeseries' + str(self.sub) + '.p'
        self.load(dat='timeseries')
        self.timeseries.file=self.path + self.name + '//Files//timeseries' + str(self.sub) + '.p'
        
        self.imagedata=Imagedata()        
        self.imagedata.file=self.path + self.name + '//Files//imagedata' + str(self.sub) + '.p'
        
        self.movie=Movie()
        self.movie.file=self.path + self.name + '/Movies/'+ str(self.sub) +'/'
        
    def forcedisp(self):
        mfile=self.infos.path+self.infos.instron_file
        data=np.loadtxt(mfile,skiprows=2,delimiter=";")
        self.load(dat='timeseries')
        self.timeseries.ti=data[:,0]
        self.timeseries.ui=data[:,1]
        self.timeseries.Fi=data[:,2]
        self.save(dat='timeseries')
    #def save(self):
    #    pickle.dump(self.infos, open(self.infos.file, "w" ) )
    #    print('------  file %s saved --------' % self.infos.file)
    def save(self,dat='infos'):
        pickle.dump(self.__getattribute__(dat), open(self.__getattribute__(dat).file, "wb" ) )
        print('------  file %s saved --------' % self.__getattribute__(dat).file)
        
    def load(self,dat='infos'):
        if os.path.isfile(self.__getattribute__(dat).file):
            #tmp = pickle.load( open(self.__getattribute__(dat).file, "r" ) )
            #if dat is "infos":
            #    for f in tmp.keys():
            #        #self.__getattribute__(dat).__setattr__(f,tmp.__getattribute__(f))
            #        self.__getattribute__(dat).__setattr__(f,tmp[f])
            #else:
            tmp = pickle.load( open(self.__getattribute__(dat).file, "rb" ) )
            for f in tmp.__dict__.keys():
                self.__getattribute__(dat).__setattr__(f,tmp.__getattribute__(f))
            print(dat+ ' reloaded')
        else:
            print('define ' + dat + ' first')
    
    def part(self,p):# trajectory
        self.Part=Part.Part()
        self.Part.p=p
        self.Part.sub=self.sub
        self.Part.name=self.name
        self.Part.DB=self.DB
        
    def time(self,t):# profile
        self.Time=Time.Time()
        self.Time.t=t
        self.Time.sub=self.sub
        self.Time.name=self.name
        self.Time.DB=self.DB
        
    def timelink(self,t,fields=["bend"],tlink="alpha"):# profile
        self.time(t)
        self.Link=Link.Link()
        req="SELECT l, sign FROM LinkProp%s WHERE (typelink='%s') ORDER BY l;" % (self.sub,tlink) 
        self.Time.fetch(req,["l","sign"])
        self.Link.links=np.copy(self.Time.data.l)
        self.Link.signs=np.copy(self.Time.data.sign)
        req="SELECT " + ",".join(fields) + " FROM Link%s WHERE (t=%s AND (" % (self.sub,t)
        for l in self.Link.links:
            req+="l=%s OR " % l
        req=req[:-4]
        req+=" )) ORDER BY l;"
        #print(req)
        self.Time.fetch(req,fields)
        
    def link(self,l,fields=["bend"],tlink="alpha"):# profile
        self.Link=Link.Link()
        self.Link.sub=self.sub
        self.Link.name=self.name
        self.Link.DB=self.DB
        req="SELECT l, sign, nx1, ny1 FROM LinkProp%s WHERE (typelink='%s') ORDER BY l;" % (self.sub,tlink) 
        self.Link.fetch(req,["l","sign","nx1","ny1"])
        self.Link.links=np.copy(self.Link.data.l)
        self.Link.signs=np.copy(self.Link.data.sign)
        self.Link.nx=np.copy(self.Link.data.nx1)
        self.Link.ny=np.copy(self.Link.data.ny1)
        req="SELECT " + ",".join(fields) + " FROM Link%s WHERE (l=%s) ORDER BY t;" % (self.sub,l)
        self.Link.fetch(req,fields)
            
    def plottrajectory(self,ax,p):
        ##### retrieve #######
        self.part(p)
        self.Part.fields(['Y','X','t'])
        i=self.Part.data.t.argsort()
        ax.plot(self.Part.data.t[i],self.Part.data.Y[i])
        ax.plot(self.Part.data.t[i],self.Part.data.X[i])
    
    def prop(self):#properties
        self.Prop=Prop.Prop()
        self.Prop.sub=self.sub
        self.Prop.name=self.name
        self.Prop.DB=self.DB
    
    def drawellipse(self,ax,p,t,color='k',linestyle='-',linewidth=1):
        ##### retrieve #######
        self.part(p)
        self.Part.fields(['a','b','phi','Y','X','t'])
        ti=np.where(self.Part.data.t==t)
        x0=self.Part.data.X[ti];y0=self.Part.data.Y[ti]
        a=self.Part.data.a[ti]/2.;b=self.Part.data.b[ti]/2.;phi=self.Part.data.phi[ti]*np.pi/180.
        Tools.plotellipse(ax,a,b,phi,x0,y0,linestyle=linestyle,linewidth=linewidth,color=color)
        
    def drawsquare(self,ax,p,t,color='k',linestyle='-',linewidth=1):
        ##### retrieve #######
        self.part(p)
        self.Part.fields(['a','b','phi','Y','X','t'])
        ti=np.where(self.Part.data.t==t)
        x0=self.Part.data.X[ti];y0=self.Part.data.Y[ti]
        a=self.infos.a;b=a;phi=self.Part.data.phi[ti]*np.pi/180.
        Tools.plotsquare(ax,a,b,phi,x0,y0,linestyle=linestyle,linewidth=linewidth,color=color)
        
    def plotsquares(self,ax,t,p=[],color="r",alpha=0.5,filter_particles=False,scaling=1,offset=[0,0]):
        self.time(t)
        self.Time.fields(['t','a','b','phi','Y','X','p'])
        self.prop()
        self.Prop.fields(["p","nx","ny"])
        #theta=(np.pi/2.-self.Time.data.phi*np.pi/180.)        
        patches = []
        if len(p) == 0:
            ps=range(len(self.Time.data.X))
        else:
            ps=np.where(self.Time.data.p == p)[0]
        for e in ps:
            if filter_particles and self.Prop.data.nx[self.Time.data.p[e]]*self.Prop.data.ny[self.Time.data.p[e]] == 0: continue
            box=ellipse2square([[self.Time.data.X[e]*scaling+offset[0],self.Time.data.Y[e]*scaling+offset[1]],[],self.Time.data.phi[e]*np.pi/180],squaresize=np.sqrt(2)/2*scaling,inte=False)
            box=np.array([list(b) for b in box])            
            sq=Polygon(box)
            patches.append(sq)
        p = PatchCollection(patches, alpha=alpha,facecolor="None",edgecolor=color,linewidth=3)
        ax.add_collection(p)
    
    def plotellipses(self,ax,t,p=[],color="r",alpha=0.5,filter_particles=False,scaling=1,offset=[0,0]):
        self.time(t)
        self.Time.fields(['t','a','b','phi','Y','X','p'])
        self.prop()
        self.Prop.fields(["p","nx","ny"])
        #theta=(np.pi/2.-self.Time.data.phi*np.pi/180.)        
        patches = []
        if len(p) == 0:
            ps=range(len(self.Time.data.X))
        else:
            ps=np.where(self.Time.data.p == p)[0]
        for e in ps:
            if filter_particles and self.Prop.data.nx[self.Time.data.p[e]]*self.Prop.data.ny[self.Time.data.p[e]] == 0: continue
            #box=ellipse2square([[self.Time.data.X[e]*scaling+offset[0],self.Time.data.Y[e]*scaling+offset[1]],[],self.Time.data.phi[e]*np.pi/180],squaresize=np.sqrt(2)/2*scaling,inte=False)
            #box=np.array([list(b) for b in box])
            el=Ellipse(xy=[self.Time.data.X[e]*scaling+offset[0],self.Time.data.Y[e]*scaling+offset[1]],
                       width=self.Time.data.a[e]*scaling,height=self.Time.data.b[e]*scaling,angle=self.Time.data.phi[e])
            #sq=Polygon(box)
            patches.append(el)
        p = PatchCollection(patches, alpha=alpha,facecolor="None",edgecolor=color,linewidth=2)
        ax.add_collection(p)
        #ax.scatter(self.Time.data.X,self.Time.data.Y) 
    def extract_bendingangles(self,tlink="beta",ax=[],tolerance1=120,tolerance2=5):
        self.load(dat='timeseries')
        self.time(0)
        self.Time.timelist()
        
        malpha=np.zeros(len(self.Time.data.timelist))
        N=np.zeros(len(self.Time.data.timelist))
        stdalpha=np.zeros(len(self.Time.data.timelist))
        N=np.zeros(len(self.Time.data.timelist))
    
        tolerance1=105
        self.link(0,tlink=tlink)
        for link,sign,nx,ny in zip(self.Link.links[:],self.Link.signs[:],self.Link.nx,self.Link.ny):
            self.link(link,fields=["t","bend"])
            
            angle1=self.Link.data.bend
            angle2=np.array([angle1[i] for i in range(len(angle1)) if isnan(angle1[i]) is False])# and np.abs(angle[i]-offset) <=np.pi/2*1.1])
            t2=np.array([self.Link.data.t[i] for i in range(len(angle1)) if isnan(angle1[i]) is False])# and np.abs(angle[i]-offset) <=np.pi/2*1.1]
            #N2=np.array([1 for i in range(len(angle1)) if isnan(angle1[i]) is False])
            angle3=jumpcorrector(np.copy(angle2),tolerance1)*sign*np.pi/180
            angle3=jumpcorrector(np.copy(angle3)*180/np.pi,tolerance1)*np.pi/180
            #            
            offset=angle3[0]#offset=np.rint(np.nanmean(angle)/np.pi)*np.pi
            if nx >=self.infos.ROI.nxmin and nx <=self.infos.ROI.nxmax and ny >=self.infos.ROI.nymin and ny<=self.infos.ROI.nymax:#True:#len(angle) == len(self.Time.data.timelist) and np.abs(angle.min()-angle.max())<=np.pi/2*1.1:
                print(len(angle3-offset),nx,ny)
                rem=np.copy(np.hstack([np.abs(np.diff(angle3))<=tolerance2*np.pi/180,True]))
                t4=t2[rem]
                angle4=angle3[rem]
                rem=np.copy(np.hstack([np.abs(np.diff(angle4))<=tolerance2*np.pi/180,True]))
                t5=t4[rem]
                angle5=angle4[rem]
                rem=np.copy(np.hstack([np.abs(np.diff(angle5))<=tolerance2*np.pi/180,True]))
                t6=t5[rem]
                angle6=angle5[rem]
                rem=np.copy(np.hstack([np.abs(np.diff(angle6))<=tolerance2*np.pi/180,True]))
                t7=t6[rem]
                angle7=angle6[rem]
                rem=np.copy(np.hstack([np.abs(np.diff(angle7))<=tolerance2*np.pi/180,True]))
                t8=t7[rem]
                angle8=angle7[rem]
                angle8=jumpcorrector(np.copy(angle8)*180/np.pi,tolerance1)*np.pi/180
                
                malpha[t8]+=angle8-offset
                stdalpha[t8]+=(angle8-offset)**2
                N[t8]+=1
                if len(ax.__dir__())>45:
                    ax.plot(t8,angle8-offset)
                
        self.timeseries.__setattr__(tlink+"_mean",malpha/N)
        self.timeseries.__setattr__(tlink+"_std",np.sqrt(stdalpha/N-malpha**2/N**2))
        self.save(dat='timeseries')
        
    def extract_bendingangles_old(self,tlink="beta",ax=[],tolerance1=120,tolerance2=5):
        self.load(dat='timeseries')
        self.time(0)
        self.Time.timelist()
        
        malpha=np.zeros(len(self.Time.data.timelist))
        N=np.zeros(len(self.Time.data.timelist))
        stdalpha=np.zeros(len(self.Time.data.timelist))
    
        self.link(0,tlink=tlink)
        for link,sign in zip(self.Link.links[:],self.Link.signs[:]):
            self.link(link)
            
            angle=-jumpcorrector(self.Link.data.bend,tolerance1,tolerance2=tolerance2)*sign*np.pi/180
            offset=angle[0]#np.rint(angle[0]/np.pi)*np.pi
            if np.abs(angle.min()-angle.max())<=np.pi/2*1.1:
                malpha+=angle-offset
                stdalpha+=(angle-offset)**2
                N+=1
                if len(ax.__dir__())>45:
                    ax.plot(angle-offset)
        self.timeseries.__setattr__(tlink+"_mean",malpha/N)
        self.timeseries.__setattr__(tlink+"_std",np.sqrt(stdalpha/N-malpha**2/N**2))
        self.save(dat='timeseries')
        
    def detect_ellipses(self,display=False,displaypic=0,save=True):
        ### Parameters of the simulation ###
        path=self.infos.path+self.infos.suffix+"/"
        pathcheck=self.infos.path+"DetectionCheck/"+self.infos.suffix+"/"
        if os.path.isdir(pathcheck) is False: os.mkdir(pathcheck)
                
        pics=os.listdir(path)        
        if True:#display is True:
            fig=plt.figure(1,figsize=(10,10))
            ax1=fig.add_subplot(111)
            #ax2=fig.add_subplot(122)
        else:
            ax1=None;ax2=None        
        num=[]
        for ff in os.listdir(path):
            if ".tif" in ff:
                num+=[int(ff[15:20])]
        num=np.array(num)
        num.sort()
        pics=num.tolist()
        
        self.load(dat='timeseries');
        self.timeseries.t=np.zeros(len(pics))
        self.timeseries.pos2=np.zeros(len(pics))
        tt=0
        self.imagedata.ellipses=[[]]*len(pics)
        if display is True:pics=[pics[displaypic]]
        #print(pics)
        for jj in pics:
            file_in= self.infos.suffix+"_%05d.tif" % (jj)
            file_out= self.infos.suffix+"_%05d.jpg" % (jj)
            #file_in= "_%05d.tif" % (jj)
            #file_out= "_%05d.jpg" % (jj)
            #file_out= "img_%05d.jpg" % (jj)            
            #Load image
            print(path + file_in)
            originalImage = cv2.imread(path + file_in)
            if originalImage is not None: #Does the image exist?
                #///////////////////// INITIAL STEP ///////////////////////
                ######## Load ####
                img = cv2.imread(path + file_in,0)#load image as grayscale image
                CLONE = np.copy(originalImage)      # to draw check picture  
                CLONE3 = np.copy(originalImage) 
                ######## Rotate and Crop ####
                x0,deltax,y0,deltay=self.infos.BBox
                imgcrop = np.copy(img[y0:y0+deltay, x0:x0+deltax])#Crop image)
                del originalImage
                #ax1.imshow(imgcrop)
                #imgray = cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
                (tmp,ThresholdedImage)= cv2.threshold(imgcrop,self.infos.ths2,255,cv2.THRESH_BINARY_INV)
                #(tmp,ThresholdedImage)= cv2.threshold(imgcrop,255-self.infos.ths2,255-self.infos.ths1,cv2.THRESH_BINARY_INV)
                del tmp
                #############################################
                ################### STEP 1 ##################
                #############################################
                ######### Erode and Dilate : open holes #####
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.erode(ThresholdedImage,kernel)
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.dilate(ThresholdedImage,kernel)
                #ax2.imshow(ThresholdedImage)
                ######## Fit to ellipses ####    
                img1=np.copy(ThresholdedImage)
                ellipses=fitellipses(img1,self.infos.perim_lim,BBox=self.infos.BBox)
                del img1
                ############# SAVE #########
                imgcrop2=np.copy(img)
                ellipses_out=save_ellipses(jj,ellipses,imgcrop2,self.infos.BBox,
                                           self.infos.minArea,self.infos.maxArea,
                                           self.infos.maxaspect,bright=255,ellipses_out=[])
                self.imagedata.ellipses[jj]=ellipses_out
                for ellipse in ellipses_out:
                    #cv2.ellipse(CLONE,ellipse,(0,255,0),2)#draw this ellips to check
                    #cv2.ellipse(CLONE,ellipse,(255,255,255),-1)
                    pt1=ellipse2point(ellipse,0,a=self.infos.square)
                    pt2=ellipse2point(ellipse,np.pi/2,a=self.infos.square)
                    pt3=ellipse2point(ellipse,np.pi,a=self.infos.square)
                    pt4=ellipse2point(ellipse,-np.pi/2,a=self.infos.square)
                    box=((pt1[0],pt1[1]), (pt2[0], pt2[1]), (pt3[0], pt3[1]), (pt4[0], pt4[1]))
                    box = np.int0(box)
                    #cv2.drawContours(CLONE,[box],0,(0,191,255),2)
                    cv2.drawContours(CLONE,[box],0,(255,255,255),-1)
                    #print("yes")
                    #cv2.ellipse(CLONE,ellipse,(125,125,125),-1)
                #if ax1 is not None:
                #    ax1.imshow(CLONE)
                #cv2.imwrite(pathcheck+file_out,CLONE)
                #############################################
                ################### STEP 2 ##################
                #############################################
                CLONE2 = cv2.cvtColor(CLONE,cv2.COLOR_RGB2GRAY)
                imgcrop = np.copy(CLONE2[y0:y0+deltay, x0:x0+deltax])#Crop image)
                #ax1.imshow(imgcrop)
                #imgray = cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
                (tmp,ThresholdedImage)= cv2.threshold(imgcrop,self.infos.ths3,255,cv2.THRESH_BINARY_INV)
                #(tmp,ThresholdedImage)= cv2.threshold(imgcrop,255-self.infos.ths2,255-self.infos.ths1,cv2.THRESH_BINARY_INV)
                del tmp
                ######### Erode and Dilate : open holes #####
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.erode(ThresholdedImage,kernel)
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.dilate(ThresholdedImage,kernel)
                #ax2.imshow(ThresholdedImage)
                ######## Fit to ellipses ####    
                img1=np.copy(ThresholdedImage)
                #ax1.imshow(ThresholdedImage)
                ellipses=fitellipses(img1,self.infos.perim_lim,BBox=self.infos.BBox)
                del img1
                ############# SAVE #########
                imgcrop2=np.copy(img)
                ellipses_out=save_ellipses(jj,ellipses,imgcrop2,self.infos.BBox,
                                           self.infos.minArea,self.infos.maxArea,
                                           self.infos.maxaspect,bright=255,ellipses_out=[])
                self.imagedata.ellipses[jj]+=ellipses_out
                for ellipse in ellipses_out:
                    #cv2.ellipse(CLONE,ellipse,(0,255,0),2)#draw this ellips to check
                    #cv2.ellipse(CLONE,ellipse,(125,125,125),-1)
                    cv2.ellipse(CLONE3,ellipse,(255,255,255),-1)
                    #print(ellipse)
                if ax1 is not None:
                    for ellipse in self.imagedata.ellipses[jj]:
                        cv2.ellipse(CLONE3,ellipse,(0,255,0),2)#draw this ellips to check
                        #cv2.ellipse(CLONE3,ellipse,(255,255,255),-1)
                    #ax1.imshow(CLONE3)
                cv2.imwrite(pathcheck+file_out,CLONE3)
                #cv2.imwrite(pathcheck+file_out,CLONE)
                #cv2.imwrite(pathcheck+file_out,ThresholdedImage)
                
                del ellipses_out,ellipses,imgcrop2,CLONE,ThresholdedImage,CLONE2,CLONE3
                    ######## Analyze #######
            tt+=1
        if save is True:self.save(dat="imagedata")
    def detect_ellipses_III(self,display=False,displaypic=0,save=True):
        ### Parameters of the simulation ###
        path=self.infos.path+self.infos.suffix+"//"
        pathcheck=self.infos.path+"DetectionCheck//"+self.infos.suffix+"//"
        if os.path.isdir(pathcheck) is False: os.mkdir(pathcheck)
                
        pics=os.listdir(path)        
        if True:#display is True:
            fig=plt.figure(1,figsize=(10,10))
            ax1=fig.add_subplot(111)
            #ax2=fig.add_subplot(122)
        else:
            ax1=None    
        num=[]
        for ff in os.listdir(path):
            if ".tif" in ff:
                num+=[int(ff[15:20])]
        num=np.array(num)
        num.sort()
        pics=num.tolist()
        
        self.load(dat='timeseries');
        self.timeseries.t=np.zeros(len(pics))
        self.timeseries.pos2=np.zeros(len(pics))
        tt=0
        self.imagedata.ellipses=[[]]*len(pics)
        if display is True:pics=[pics[displaypic]]
        #print(pics)
        for jj in pics:
            file_in= self.infos.suffix+"_%05d.tif" % (jj)
            file_out= self.infos.suffix+"_%05d.jpg" % (jj)
            file_out2= self.infos.suffix+"check_%05d.jpg" % (jj)
            #Load image
            print(path + file_in)
            originalImage = cv2.imread(path + file_in)
            if originalImage is not None: #Does the image exist?
                #///////////////////// INITIAL STEP ///////////////////////
                ######## Load ####
                img = cv2.imread(path + file_in,0)#load image as grayscale image
                CLONE = np.copy(originalImage)      # to draw check picture  
                CLONE_OUT = np.copy(originalImage) 
                ######## Rotate and Crop ####
                x0,deltax,y0,deltay=self.infos.BBox
                imgcrop = np.copy(img[y0:y0+deltay, x0:x0+deltax])#Crop image)
                del originalImage
                (tmp,ThresholdedImage)= cv2.threshold(imgcrop,self.infos.ths2,255,cv2.THRESH_BINARY_INV)
                del tmp
                #############################################
                ################### STEP 1 ##################
                #############################################
                ######### Erode and Dilate : open holes #####
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.erode(ThresholdedImage,kernel)
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.dilate(ThresholdedImage,kernel)
                ######## Fit to ellipses ####    
                img1=np.copy(ThresholdedImage)
                ellipses=fitellipses(img1,self.infos.perim_lim,BBox=self.infos.BBox)
                del img1
                ############# SAVE #########
                imgcrop2=np.copy(img)
                ellipses_out=save_ellipses(jj,ellipses,imgcrop2,self.infos.BBox,
                                           self.infos.minArea,self.infos.maxArea,
                                           self.infos.maxaspect,bright=255,ellipses_out=[])
                self.imagedata.ellipses[jj]=ellipses_out
                ############# OVERLAY SQUARES (for next steps) AND ELLIPSES (for check) #########                
                for ellipse in ellipses_out:
                    box=ellipse2square(ellipse,squaresize=self.infos.square)
                    cv2.drawContours(CLONE,[box],0,(255,255,255),-1)
                if ax1 is not None:
                    for ellipse in ellipses_out:
                        cv2.ellipse(CLONE_OUT,ellipse,(255,125,0),2)
                #cv2.imwrite(pathcheck+file_out,CLONE)
                #############################################
                ################### STEP 2 ##################
                ############################################# 
                #CLONE2 = cv2.cvtColor(np.copy(CLONE),cv2.COLOR_RGB2GRAY) 
                #print(CLONE.shape)
                CLONE2 = cv2.cvtColor(CLONE,cv2.COLOR_RGB2GRAY) 
                imgcrop = np.copy(CLONE2[y0:y0+deltay, x0:x0+deltax])#Crop image)
                (tmp,ThresholdedImage)= cv2.threshold(imgcrop,self.infos.ths3,255,cv2.THRESH_BINARY_INV)
                del tmp                
                ######### Erode and Dilate : open holes #####
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.erode(ThresholdedImage,kernel)
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.dilate(ThresholdedImage,kernel)
                ######## Fit to ellipses ####    
                img1=np.copy(ThresholdedImage)
                ellipses=fitellipses(img1,self.infos.perim_lim,BBox=self.infos.BBox)
                del img1
                ############# SAVE #########
                imgcrop2=np.copy(img)
                ellipses_out=save_ellipses(jj,ellipses,imgcrop2,self.infos.BBox,
                                           self.infos.minArea,self.infos.maxArea,
                                           self.infos.maxaspect,bright=255,ellipses_out=[])
                self.imagedata.ellipses[jj]+=ellipses_out
                ############# OVERLAY SQUARES (for next steps) AND ELLIPSES (for check) #########
                for ellipse in ellipses_out:
                    box=ellipse2square(ellipse,squaresize=self.infos.square)
                    cv2.drawContours(CLONE2,[box],0,(255,255,255),-1)
                if ax1 is not None:
                    for ellipse in ellipses_out:
                        cv2.ellipse(CLONE_OUT,ellipse,(0,255,0),2)
                #print(CLONE2.shape)
                
                #############################################
                ################### STEP 3 ##################
                ############################################# 
                CLONE3 = np.copy(CLONE2)#cv2.cvtColor(np.copy(CLONE2),cv2.COLOR_RGB2GRAY)
                #print(CLONE3.dtype)
                #CLONE = cv2.cvtColor(np.copy(CLONE),cv2.COLOR_RGB2GRAY)
                #CLONE = cv2.cvtColor(CLONE,cv2.COLOR_RGB2GRAY) 
                imgcrop = np.copy(CLONE3[y0:y0+deltay, x0:x0+deltax])#Crop image)
                (tmp,ThresholdedImage)= cv2.threshold(imgcrop,self.infos.ths1,255,cv2.THRESH_BINARY_INV)
                del tmp       
                
                ######### Erode and Dilate : open holes #####
                kernel = np.ones((15,15),np.uint8)
                ThresholdedImage=cv2.dilate(ThresholdedImage,kernel)
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.erode(ThresholdedImage,kernel)
                #cv2.imwrite(pathcheck+file_out,ThresholdedImage)
                ######## Fit to ellipses ####    
                img1=np.copy(ThresholdedImage)
                ellipses=fitellipses(img1,self.infos.perim_lim,BBox=self.infos.BBox)
                del img1
                ############# SAVE #########
                imgcrop2=np.copy(img)
                ellipses_out=save_ellipses(jj,ellipses,imgcrop2,self.infos.BBox,
                                           self.infos.minArea,self.infos.maxArea,
                                           self.infos.maxaspect,bright=255,ellipses_out=[])
                self.imagedata.ellipses[jj]+=ellipses_out
                if ax1 is not None:
                    for ellipse in ellipses_out:
                        cv2.ellipse(CLONE_OUT,ellipse,(0,125,255),2)
                #print(CLONE3.shape)
                ############# DETECTION CHECK #########
                cv2.imwrite(pathcheck+file_out2,CLONE2)
                cv2.imwrite(pathcheck+file_out,CLONE_OUT)
                #cv2.imwrite(pathcheck+file_out,CLONE2)
                #cv2.imwrite(pathcheck+file_out,ThresholdedImage)
                
                del ellipses_out,ellipses,imgcrop2,CLONE,CLONE2,CLONE3,ThresholdedImage,CLONE_OUT
                    ######## Analyze #######
            tt+=1
        if save is True:self.save(dat="imagedata")
    def track_ellipses(self):
        self.load(dat="imagedata")
        #tpos=np.zeros(7)
        pathfile=self.path + self.name + '/Files/imagedata_raw' + str(self.sub) + '.txt'
        text_file = open(pathfile, 'wb')
        #Header of the text file
        text_file.writelines([b't\t p\t x0\t y0\t a\t b\t angle\n'])
        #start = time.time()
        for jj in range(len(self.imagedata.ellipses)):
            for kk in range(len(self.imagedata.ellipses[jj])):
                text_file.write(b"%i\t%i\t%f\t%f\t%f\t%f\t%f\n" % (jj,kk,self.imagedata.ellipses[jj][kk][0][0],self.imagedata.ellipses[jj][kk][0][1],self.imagedata.ellipses[jj][kk][1][0],self.imagedata.ellipses[jj][kk][1][1],self.imagedata.ellipses[jj][kk][2]))
                #tpos= np.vstack((tpos,np.array([jj,kk,self.imagedata.ellipses[jj][kk][0][0],self.imagedata.ellipses[jj][kk][0][1],self.imagedata.ellipses[jj][kk][1][0],self.imagedata.ellipses[jj][kk][1][1],self.imagedata.ellipses[jj][kk][2]])))     
            #if np.mod(jj,100) == 0: 
            #    end = time.time()
            #    print("timestep %i / %i, %i sec" % (jj,len(self.imagedata.ellipses),end-start))
            #    start = time.time()
        text_file.close()
        
        self.imagedata.tpos=np.loadtxt(pathfile,skiprows=1)
        ths=self.infos.square/2.
        print("Threshold: %s" % ths)
        self.imagedata.P=tracking(self.imagedata.tpos,usesparse=True,display=True,mem=2,delta=ths)
        
        self.imagedata.T=[[]]*len(self.imagedata.ellipses)
        for t in range(len(self.imagedata.ellipses)):
            tpos=np.zeros(6)
            for p in range(len(self.imagedata.P)):
                if np.size(self.imagedata.P[p][0])>1:
                    ix=np.where(t== self.imagedata.P[p][:,0])
                    if len(ix[0])>0:
                        tmp=np.hstack((p,self.imagedata.P[p][ix[0],1:][0]))
                        tpos=np.vstack((tpos,tmp))
            if np.size(tpos[0])>1: self.imagedata.T[t]=tpos[1:,:]
        self.save(dat="imagedata")
        
    def correct_displacement(self):
        self.load(dat='timeseries')
        x=np.linspace(0,2048,200)
        f1=interp1d(np.polyval(self.infos.distortion,x),x)
        self.timeseries.pos2_corrected=f1(self.infos.BBox[0]+self.timeseries.pos2)*self.infos.pix2mm
        self.save(dat='timeseries')
        

    def make_video2(self,displaypic=0):
        if os.path.isdir(self.movie.file) is False: os.mkdir(self.movie.file)
        if os.path.isdir(self.movie.file) is False: os.mkdir(self.movie.file)
        ### Parameters of the simulation ###
        path=self.path + self.name +"/"+self.infos.suffix1+"//"+self.infos.suffix2+"//"
        pics=os.listdir(path)        

        num=[]
        for ff in os.listdir(path):
            if ".png" in ff:
                num+=[int(ff[15:20])]
        num=np.array(num)
        num.sort()
        pics=num.tolist()
        pics=[pics[displaypic]]
        #self.load(dat='timeseries');
        
        jj0=int(np.round(self.infos.Fzero_index))
        file_in= self.infos.suffix2+"_%05d.png" % (jj0)
        originalImage = cv2.imread(path + file_in,0)
        REFERENCE = np.copy(originalImage)
        for jj in pics:
            file_in= self.infos.suffix2+"_%05d.png" % (jj)
            print(path + file_in)
            originalImage = cv2.imread(path + file_in,0)
            if originalImage is not None: #Does the image exist?
                #///////////////////// INITIAL STEP ///////////////////////
                ######## Load ####
                img = cv2.imread(path + file_in)
                CLONE = np.copy(originalImage)      # to draw check picture  
                self.movie.img=img
                self.movie.img_diff=CLONE/2+128-REFERENCE/2

#Define custom links depending on the positions of the squares and the rank of the hierarchical structure.
def structure_rankII(nx,ny,nx2,ny2):
    sign=0
    if (np.mod(nx,3)==2 and np.mod(ny,3)==2) :
        typelink="beta";sign=1
    elif (np.mod(nx2,3)==2 and np.mod(ny2,3)==2):
        typelink="beta";sign=-1   
    else:
        typelink="alpha";
        sign=2*np.mod(np.floor(nx/3)+np.floor(ny/3),2)-1
    return sign, typelink
    
def structure_rankII_chirality(nx,ny,nx2,ny2):
    sign=0
    if (np.mod(nx,3)==2 and np.mod(ny,3)==2) :
        typelink="beta";sign=(2*np.mod(np.floor(nx/3)+np.floor(ny/3),2)-1)
    elif (np.mod(nx2,3)==2 and np.mod(ny2,3)==2):
        typelink="beta";sign=-(2*np.mod(np.floor(nx2/3)+np.floor(ny2/3),2)-1)   
    else:
        typelink="alpha";
        sign=2*np.mod(np.floor(nx/3)+np.floor(ny/3),2)-1
    return sign, typelink
    
def structure_rankIII(nx,ny,nx2,ny2):
    sign=0    
    typelink="TBD"
    if (np.mod(nx,3)==2 and np.mod(ny,3)==2):
        typelink="gamma";sign=-1
    elif (np.mod(nx2,3)==2 and np.mod(ny2,3)==2):
        typelink="gamma";sign=1
    elif (np.mod(nx,3)==0 and np.mod(ny,3)==2) or (np.mod(nx,3)==2 and np.mod(ny,3)==0):
        typelink="beta";
        if (np.mod(nx,9)==6 or np.mod(ny,9)==6):
            sign=-1
        else:
            sign=1
    if (np.mod(nx,9)==0 and np.mod(ny,9)==5) or (np.mod(nx,9)==5 and np.mod(ny,9)==0):
        typelink="alpha";
        sign=-(2*np.mod(np.floor(nx/9)+np.floor(ny/9),2)-1)
    #else:
    #    typelink="beta";
    #    sign=2*np.mod(np.floor(nx/3)+np.floor(ny/3),2)-1
        
    return sign, typelink

structure_hierarchical={"20170214b": structure_rankII,
                        "20170214": structure_rankIII}

class CreateDB(Focus):
    
    def __init__(self,name,sub):
        self.pref='/'
        self.name=name
        self.sub=sub
        self.DB=Database.DB()
        self.DB.connect(self.name)
        self.tmp=Tmp()
        
        self.infos=Infos()
        self.infos.file=self.path + self.name + '/Files/infos' + str(self.sub) + '.p'
        
        self.imagedata=Imagedata()        
        self.imagedata.file=self.path + self.name + '/Files/imagedata' + str(self.sub) + '.p'
        
        self.timeseries=Timeseries()
        self.timeseries.file=self.path + self.name + '/Files/timeseries' + str(self.sub) + '.p'
        
    def fetch(self,req,Fields):#
         #print '\r' + req
         self.DB.cursor.execute(req)
         ps=[]
         [ps.append(p) for p in self.DB.cursor]
         for ind,f in zip(range(len(Fields)),Fields):#f in Fields:
            #self.tmp.__setattr__(f,np.array([line[f] for line in self.DB.cursor]))
            self.tmp.__setattr__(f,np.array([p[ind] for p in ps]))
    
    def calibration(self):       
        try:
            mfile=self.path + self.name+ self.pref + str(self.sub)+'/calibration.txt'
            print(mfile)
            data=np.loadtxt(mfile,skiprows=1)
            self.infos.pix2mm=data[0]
        except: 
            traceback.print_exc()
            self.infos.pix2mm=1.
            
    def writetxtDATAPROP(self):
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)
        
        fid=open(tmpdir + '/Data.txt','w')
        fid2=open(tmpdir + '/Prop.txt','w')
        for p,pp in zip(range(len(self.imagedata.P)),self.imagedata.P):
            try:
                t=pp[:,0];X=pp[:,1]*self.infos.pix2mm/4.5;Y=pp[:,2]*self.infos.pix2mm/4.5;
                a=pp[:,3]*self.infos.pix2mm/4.5;b=pp[:,4]*self.infos.pix2mm/4.5;phi=pp[:,5];
                tt=0
                for i in range(len(t)):
                    fid.write('%d;%d;%f;%f;%f;%f;%f\n'% (p,t[i],X[i],Y[i],a[i],b[i],phi[i]))
                    fid2.write('%d\n'% (p))
                    tt+=1
            except:
                print(p)
                t=pp[0];X=pp[1]*self.infos.pix2mm/4.5;Y=pp[2]*self.infos.pix2mm/4.5;
                a=pp[3]*self.infos.pix2mm/4.5;b=pp[4]*self.infos.pix2mm/4.5;phi=pp[5];
                fid.write('%d;%d;%f;%f;%f;%f;%f\n'% (p,t,X,Y,a,b,phi))
                fid2.write('%d\n'% (p))
        fid.close()
        fid2.close()
        
    def writetxtLINK(self):
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Link.txt','w')
        fid2=open(tmpdir + '/LinkProp.txt','w')

        #all the possible locations in the grid        
        req='SELECT DISTINCT nx FROM Prop%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx'])
        nxs=np.copy(self.tmp.nx)
        req='SELECT DISTINCT ny FROM Prop%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.copy(self.tmp.ny) 
        #build links
        link=0
        for nx in nxs[:-1]:
            for ny in nys[:-1]:
                req='SELECT p FROM Prop%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
                #print(req)
                self.fetch(req,['p'])
                if len(self.tmp.p) != 1: continue
                p1=np.copy(self.tmp.p)[0]
                
                req='SELECT p,nx,ny FROM Prop%s WHERE ((nx=%s AND ny=%s) OR (nx=%s AND ny=%s));' % (self.sub,nx+1,ny,nx,ny+1)
                #print(req)
                self.fetch(req,['p','nx','ny'])
                ps=np.copy(self.tmp.p)
                nxs2=np.copy(self.tmp.nx);nys2=np.copy(self.tmp.ny)
                if len(ps) == 0: continue
                
                self.part(p1)
                self.Part.fields(["X","Y","t","phi"])
                Xp1=np.copy(self.Part.data.X);Yp1=np.copy(self.Part.data.Y);tp1=np.copy(self.Part.data.t);phip1=np.copy(self.Part.data.phi)

                for p2,nx2,ny2 in zip(ps,nxs2,nys2):                
                    self.part(p2)
                    self.Part.fields(["X","Y","t","phi"])
                    Xp2=np.copy(self.Part.data.X);Yp2=np.copy(self.Part.data.Y);tp2=np.copy(self.Part.data.t);phip2=np.copy(self.Part.data.phi)
                    ts=list(set(tp1).intersection(tp2))
                    
                    sign,typelink=structure_hierarchical[self.name](nx,ny,nx2,ny2)
                    if self.name == "20170214b" and self.sub == 9: sign,typelink=structure_rankII_chirality(nx,ny,nx2,ny2)
                    #if self.name == "20170214" and self.sub == 9: sign,typelink=structure_rankII_chirality(nx,ny,nx2,ny2)
                    fid2.write('%d;%d;%d;%d;%d;%d;%d;%s;%d\n'% (link,p1,p2,nx,ny,nx2,ny2,typelink,sign))
                    for tt in ts:
                        ix1=np.where(tp1== tt)
                        ix2=np.where(tp2== tt)
                        distance=np.sqrt((Xp1[ix1]-Xp2[ix2])**2+(Yp1[ix1]-Yp2[ix2])**2)
                        angle=np.arctan2(Xp1[ix1]-Xp2[ix2],Yp1[ix1]-Yp2[ix2])
                        bending=phip2[ix2]-phip1[ix1]  
                        fid.write('%d;%d;%f;%f;%f\n'% (link,tt,distance,angle,bending))
                    link+=1
        fid.close()
        fid2.close()

        req='DROP TABLE Link' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE LinkProp' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass        
        
        req='CREATE TABLE Link' + str(self.sub) + ' ( l int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , delta float NOT  NULL , ori float NOT NULL , bend float NOT NULL , PRIMARY KEY ( t ,  l )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE LinkProp' + str(self.sub) + ' ( l int( 8  )  NOT  NULL , p1 int( 8  )  NOT  NULL , p2 int( 8  )  NOT  NULL ,  nx1 int( 8  )  NOT  NULL , ny1 int( 8  )  NOT  NULL , nx2 int( 8  )  NOT  NULL , ny2 int( 8  )  NOT  NULL ,typelink char( 5 ) NOT  NULL, sign int( 2  )  NOT  NULL , PRIMARY KEY ( l , p1 , p2 ,typelink )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Link.txt' INTO TABLE Link" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/LinkProp.txt' INTO TABLE LinkProp" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        
    def createTables(self):
        req='CREATE TABLE Data' + str(self.sub) + ' ( p int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL , a float NOT NULL , b float NOT NULL , phi float NOT NULL , PRIMARY KEY (  t ,  p  )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE Prop' + str(self.sub) + ' ( p int( 8  )  NOT  NULL , nx int( 8  ), ny int( 8  )  NOT  NULL , PRIMARY KEY (  p, nx, ny )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)

    def dropTables(self):
        req='DROP TABLE Prop' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE Data' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass
    def pushtoDB(self):
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Data.txt' INTO TABLE Data" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Prop.txt' INTO TABLE Prop" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        
    def makeDB(self):
        self.writetxtDATAPROP()
        self.dropTables()
        self.createTables()
        self.pushtoDB()


### class DBReshape with inheritence Sheet, to modify the Database ###
class DBReshape(CreateDB):

    def __init__(self,name,sub):
        self.name=name
        self.sub=sub
        self.DB=Database.DB()
        self.DB.connect(self.name)
        self.tmp=Tmp()
        
        self.infos=Infos()
        self.infos.file=self.path + self.name + '/Files/infos' + str(self.sub) + '.p'

    def compress_indexes(self):
        ############ compress the indices such that the index corresponds to p #########
        req='SELECT DISTINCT p FROM Data' + str(self.sub) + ' ORDER BY p;'
        self.fetch(req,['p'])
        p0=np.copy(self.tmp.p)
        print("re-indexing of the table")
        for i in range(len(p0)):
            t=0
            if p0[i] != i:
                req='UPDATE Data' + str(self.sub) + ' SET p=' + str(i) + ' WHERE p=' + str(p0[i]) + ';'
                self.DB.cursor.execute(req)
                #print(req)
                req='UPDATE Prop' + str(self.sub) + ' SET p=' + str(i) + ' WHERE p=' + str(p0[i]) + ';'
                #print(req)
                self.DB.cursor.execute(req)#print(req)
                t+=1
        if t==0:print("nothing has been re-indexed")
    
    def group(self,ps):
        ############ Merge the p values of several particles #########
        print("Check if not overlapping time steps")
        req='SELECT p,t FROM Data' + str(self.sub) + ' WHERE ('
        for p in ps:
            req+='p=' + str(p) + ' OR '
        req=req[:-4:]
        req+=');'
        #print(req)
        self.fetch(req,['p','t'])
        if len(self.tmp.t)!=len(set(self.tmp.t)):
            print("overlapping time steps: Aborting")
        else:
            print("Grouping...")
            ### changing Data ###
            pnew=ps.min()
            req='UPDATE Data' + str(self.sub) + ' SET p=' + str(pnew) + ' WHERE ('
            for p in ps:
                req+='p=' + str(p) + ' OR '
            req=req[:-4:]
            req+=');'
            print(req)
            self.DB.cursor.execute(req)
            ### changing Prop ###
            pdel=ps[np.where(ps != pnew)]
            if len(pdel)!=0:
                req='DELETE From Prop' + str(self.sub) + ' WHERE ('
                for p in pdel:
                    req+='p=' + str(p) + ' OR '
                req=req[:-4:]
                req+=');'
                print(req)
                self.DB.cursor.execute(req)
                
    def group2(self,ps):
        ############ Merge the p values of several particles #########
        print("Check if not overlapping time steps")
        req='SELECT p,t FROM Data%s WHERE p=%s;' % (self.sub,ps[0])
        print(req)
        self.fetch(req,['p','t'])
        t1=np.copy(self.tmp.t)
        req='SELECT p,t FROM Data%s WHERE p=%s;' % (self.sub,ps[1])
        print(req)
        self.fetch(req,['p','t'])
        t2=np.copy(self.tmp.t)
        overlap=set(t1).intersection(t2)
        if len(overlap)>0:
            print("overlapping time step at times %s\n" % list(overlap))
            rep=input("Force overlap? particle label with be moved from %s to %s except at the overlapping times (y/n)\n" % (ps[1],ps[0]))
            if rep == "y":
                req='UPDATE Data%s SET p=%s WHERE (p=%s AND (' % (self.sub,ps[0],ps[1])
                for t in list(set(t2)-overlap):
                    req+='t=' + str(t) + ' OR '
                req=req[:-4:]
                req+='));'
                print(req)
                self.DB.cursor.execute(req)
                self.compress_indexes()
        else:
            print("Grouping...")
            ### changing Data ###
            pnew=ps.min()
            req='UPDATE Data' + str(self.sub) + ' SET p=' + str(pnew) + ' WHERE ('
            for p in ps:
                req+='p=' + str(p) + ' OR '
            req=req[:-4:]
            req+=');'
            print(req)
            self.DB.cursor.execute(req)
            ### changing Prop ###
            #pdel=ps[find(ps != pnew)]
            pdel=ps[np.where(ps != pnew)]
            if len(pdel)!=0:
                req='DELETE From Prop' + str(self.sub) + ' WHERE ('
                for p in pdel:
                    req+='p=' + str(p) + ' OR '
                req=req[:-4:]
                req+=');'
                print(req)
                self.DB.cursor.execute(req)

            
    def split(self,p,t):
        ############ split two trajectories and assign two distinct values of p #########
        req='SELECT DISTINCT p FROM Data' + str(self.sub) + ' ORDER BY p;'
        self.fetch(req,['p'])
        p0=np.copy(self.tmp.p)
        print("Splitting trajectories into two distinct ones at time t %d" % t)
        pnew=p0.max()+1
        #req='UPDATE Data' + str(self.sub) + ' SET p=' + str(pnew) + ' WHERE (p=' +  str(p) + ' AND t>=' + str(t) +');'
        #req="UPDATE Data%(sub)s SET p=%(pnew)s WHERE (p=%(p)s AND t>=%(t)s);"
        req="UPDATE Data%s SET p=%s WHERE (p=%s AND t>=%s);" % (self.sub,pnew,p,t)
        print(req)
        self.DB.cursor.execute(req)
        req='INSERT INTO Prop' + str(self.sub) + ' (p,nx,ny) VALUES(' + str(pnew) + ',0,0)'
        print(req)
        self.DB.cursor.execute(req)

    def remove(self,p):
        sub=self.sub
        ############ split two trajectories and assign two distinct values of p #########
        print("Erasing particle %(p)s from databases")
        #req='DELETE FROM Data%(sub)s WHERE p=%(p)s;'
        req='DELETE FROM Data%s WHERE p=%s;' % (sub,p)
        #print(req % locals())
        print(req)
        #self.DB.cursor.execute(req, locals())
        self.DB.cursor.execute(req)
        #req='DELETE FROM Prop%(sub)s WHERE p=%(p)s;'
        req='DELETE FROM Prop%s WHERE p=%s;' % (sub,p)
        #print(req % locals())
        #self.DB.cursor.execute(req, locals())
        print(req)
        self.DB.cursor.execute(req)
        
    def remove_batch(self,pdel):
        #pdel=ps[np.where(ps != pnew)]
        if len(pdel)!=0:
            req='DELETE From Prop' + str(self.sub) + ' WHERE ('
            for p in pdel:
                req+='p=' + str(p) + ' OR '
            req=req[:-4:]
            req+=');'
            print(req)
            self.DB.cursor.execute(req)
            req='DELETE From Data' + str(self.sub) + ' WHERE ('
            for p in pdel:
                req+='p=' + str(p) + ' OR '
            req=req[:-4:]
            req+=');'
            print(req)
            self.DB.cursor.execute(req)

    def sortradius(self):
        sub=self.sub
        ############ compress the indices such that the index corresponds to p #########
        self.prop()
        self.Prop.fields(['p','R'])
        p0=np.copy(self.Prop.data.p)
        R0=np.copy(self.Prop.data.R)
        print("re-indexing of the table")
        for i in range(len(p0)):
            p=p0[i]
            if R0[i]>7: r=2
            else: r=1
            req='UPDATE Prop%(sub)s SET R=%(r)s WHERE p=%(p)s;'
            print(req % locals())
            self.DB.cursor.execute(req,locals())

    def swapradius(self,p):
        sub=self.sub
        ############ compress the indices such that the index corresponds to p #########
        self.prop()
        self.Prop.fields(['p','R'])
        ix=np.where(self.Prop.data.p==p)
        r=self.Prop.data.R[ix]
        if self.Prop.data.R[ix]==1:   r=2
        elif self.Prop.data.R[ix]==2: r=1
        req='UPDATE Prop%(sub)s SET R=%(r)s WHERE p=%(p)s;'
        print(req % locals())
        self.DB.cursor.execute(req,locals())
    
    def updatenx(self,pu,nx):
        sub=self.sub
        ############ compress the indices such that the index corresponds to p #########
        self.prop()
        self.Prop.fields(['p'])
        req='UPDATE Prop%s SET nx=%s  WHERE (' % (sub,nx)
        for p in list(pu):
            req+='p=' + str(p) + ' OR '
        req=req[:-4:]
        req+=');'
        print(req)
        #print(req % locals())
        self.DB.cursor.execute(req)

    def updateny(self,pu,ny):
        sub=self.sub
        ############ compress the indices such that the index corresponds to p #########
        self.prop()
        self.Prop.fields(['p'])
        req='UPDATE Prop%s SET ny=%s  WHERE (' % (sub,ny)
        for p in list(pu):
            req+='p=' + str(p) + ' OR '
        req=req[:-4:]
        req+=');'
        print(req)
        #print(req % locals())
        self.DB.cursor.execute(req)

    def remove_particle(self,ax1):
        X1=self.Time.data.X;Y1=self.Time.data.Y;p1=self.Time.data.p;t=self.Time.data.t[0]
        pi1,ix1=Tools.pick(p1,X1,Y1)
        self.plotsquares(ax1,t,p=[pi1],color="b")
        #F.drawellipse(ax1,pi1,t1[ix1],color='b',linestyle='-',linewidth=1)
        x=plt.ginput(1)
        conf=input('continue? (1=y/0=n)\n')
        if conf  == "1":
            self.remove(pi1)
            self.compress_indexes()
    
    def remove_batch_particles(self,ax1):
        X1=self.Time.data.X;Y1=self.Time.data.Y;p1=self.Time.data.p;t=self.Time.data.t[0]
        go = "1"
        ps=[]
        while go == "1":
            pi1,ix1=Tools.pick(p1,X1,Y1)
            ps+=[pi1]
            self.plotsquares(ax1,t,p=[pi1],color="b")
            go= input("keep on selecting particles? (1=y/0=n)\n")
        x=plt.ginput(1)
        conf=input('This will erase the particles, continue? (1=y/0=n)\n')
        if conf  == "1":
            self.remove_batch(ps)
            self.compress_indexes()
    

    def plot_trajectories(self,ax1):
        X1=np.copy(self.Time.data.X);Y1=np.copy(self.Time.data.Y);p1=np.copy(self.Time.data.p);t=np.copy(self.Time.data.t[0])
        
        tb=int(input("time frame 2?")) 
        self.time(tb)
        self.Time.timelist()
        self.plotsquares(ax1,tb,color="b")
        X2=np.copy(self.Time.data.X);Y2=np.copy(self.Time.data.Y);p2=np.copy(self.Time.data.p)#;t2=F.Time.data.t
        
        
        print("test")
        for pp,id1,xx1,yy1 in zip(p1,range(len(p1)),X1,Y1):
            id2=np.where(p2 == pp)[0]
            if len(id2) ==1:
                xx2=X2[id2];yy2=Y2[id2]
                print([xx1,xx2,yy1,yy2])
                ax1.plot([xx1,xx2],[yy1,yy2],color="k",lw=1)        
            elif len(id2) > 1:
                print("Strange, particle %s found twice at time step %s" % (pp,tb))
        
        pi1,ix1=Tools.pick(p1,X1,Y1)
        self.plotsquares(ax1,t,p=[pi1],color="r",alpha=1)
        pi2,ix2=Tools.pick(p2,X2,Y2)
        self.plotsquares(ax1,tb,p=[pi2],color="b",alpha=1)
        
        self.part(pi1)
        self.Part.fields(["X","Y","t"])
        Xp1=np.copy(self.Part.data.X);Yp1=np.copy(self.Part.data.Y);tp1=np.copy(self.Part.data.t)
        self.part(pi2)
        self.Part.fields(["X","Y","t"])
        Xp2=np.copy(self.Part.data.X);Yp2=np.copy(self.Part.data.Y);tp2=np.copy(self.Part.data.t)
        #ax1.plot(Xp1,Yp1,color="r")
        #ax1.plot(Xp2,Yp2,color="b")
        
        self.fig2=plt.figure(2,figsize=(8,8))
        self.ax2=self.fig2.add_axes([0.1,0.1,.8,.4])
        self.ax2b=self.fig2.add_axes([0.1,0.55,.8,.4])
        
        self.ax2.plot(tp1,Xp1,color="r")        
        self.ax2.plot(tp2,Xp2,color="b")
        self.ax2b.plot(tp1,Yp1,color="r")        
        self.ax2b.plot(tp2,Yp2,color="b")   
        
        self.pi1=pi1
        self.pi2=pi2
        x=plt.ginput(1)
    def group_particles(self,ax1):
        self.plot_trajectories(ax1)
        conf=input('this will group the two trajectories of particles %s and %s, continue? (1=y/0=n)' % (self.pi1,self.pi2))
        if conf=="1":
            self.group2(np.array([self.pi1,self.pi2]))
        self.fig2.clear()
        self.__delattr__("pi1");self.__delattr__("pi2")
        
    def split_trajectory(self,ax1):
        self.plot_trajectories(ax1)
        x=plt.ginput(1)
        tsplit=int(np.floor(x[0][0]))
        print(tsplit)
        self.ax2.axvline(x=tsplit)
        self.ax2b.axvline(x=tsplit)
        x1=plt.ginput(1)
        rep=input('this will split the trajectory of particle %s, continue? (1=y/0=n)' % self.pi1)
        if rep==1:# and self.pi1==self.pi2: 
            self.split(self.pi1,tsplit)
        #self.fig2.clear()
        #self.__delattr__("pi1");self.__delattr__("pi2")
        
    def assign_ny(self,ax1):
        X1=self.Time.data.X;Y1=self.Time.data.Y;p1=self.Time.data.p;
        pi1,ix1=Tools.pickrectangle(ax1,p1,X1,Y1)
        
        x=plt.ginput()
        conf=int(input('which ny?\n'))
        self.updateny(pi1,conf)  
    
    def assign_nx(self,ax1):
        X1=self.Time.data.X;Y1=self.Time.data.Y;p1=self.Time.data.p;
        pi1,ix1=Tools.pickrectangle(ax1,p1,X1,Y1)
        
        x=plt.ginput()
        conf=int(input('which nx?\n'))
        self.updatenx(pi1,conf)
        
    def assign_nxny_auto(self,ax1):
        X1=self.Time.data.X;Y1=self.Time.data.Y;p1=self.Time.data.p;
        
        N=27
        for i in range(N):
            x1=(i+1)*np.sqrt(2)
            x2=(i+2)*np.sqrt(2)
            ax1.axvline(x=x1)
            ax1.axhline(y=x1)
            
            condX1=np.double((X1>=x1)==True)
            condX2=np.double((X1<=x2)==True)
                                    
            ixX=find(condX1+condX2==2)
            ixX=np.array(ixX,int)
            poutX=p1[ixX]        
            self.updatenx(poutX,i+1)
            
            condY1=np.double((Y1>=x1)==True)
            condY2=np.double((Y1<=x2)==True)
            
            ixY=find(condY1+condY2==2)
            ixY=np.array(ixY,int)
            poutY=p1[ixY]        
            self.updateny(poutY,i+1)
        x=plt.ginput()
