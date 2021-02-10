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
#import imutils
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
from scipy.signal import savgol_filter
from sys import platform as _platform
from matplotlib.mlab import find
from math import isnan
import copy
import pandas as pd

import numpy as np
from sys import platform as _platform



path_sys={"darwin":"/Users/coulais/science/David_Corentin/Multimode/",
          "win32":r"C:/Users/David/Documents/0_PhD/15_Viscoelastic_Metamaterial/Multimode2_v02/",
          "linux":"/home/aleksi/Desktop/Metacombinatorial"}

sys.path.append(path_sys[_platform]+'PyLab')
from GraphicTools import fitellipses,save_ellipses,plot_ellipses, tracking,ellipse2square,ellipse2point
from Data import Database,Part,Time,Prop,Link
from Graphics import Tools

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

### big class Focus ###
class Focus(object):
    path=path_sys[_platform]
    print(path)
    
    name=''
    sub=np.NaN
    def __init__(self,name,sub,display=False):
        self.path=path_sys[_platform]+'Experiments/'#"/Volumes/data/AMOLF/groups/hecke-group/Corentin/rawdata/"        
        self.name=name
        self.sub=sub      
               
        self.DB=Database.DB()#to comment if no mysql server 
        self.DB.connect(self.name)#to comment if no mysql server 
        self.tmp=Tmp()
        
        self.infos=Infos()
        self.infos.file=self.path + self.name + '/Files/infos' + str(self.sub) + '.p'
        self.load(dat='infos')
        self.infos.file=self.path + self.name + '/Files/infos' + str(self.sub) + '.p'
        #self.infos.ROI=ROI()
        
        self.timeseries=Timeseries()
        self.timeseries.file=self.path + self.name + '/Files/timeseries' + str(self.sub) + '.p'
        self.load(dat='timeseries')
        self.timeseries.file=self.path + self.name + '/Files/timeseries' + str(self.sub) + '.p'
        
        self.imagedata=Imagedata()        
        self.imagedata.file=self.path + self.name + '/Files/imagedata' + str(self.sub) + '.p'
        
        self.movie=Movie()
        self.movie.file=self.path + self.name + '/Movies/'+ str(self.sub) +'/'
        
    def forcedisp(self):
        mfile=self.infos.path_rawdata+self.infos.instron_file
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
    
    def prop(self):#properties
        self.Prop=Prop.Prop()
        self.Prop.sub=self.sub
        self.Prop.name=self.name
        self.Prop.DB=self.DB
    
    def fetch(self,req,Fields):#
         #print '\r' + req
         self.DB.cursor.execute(req)
         ps=[]
         [ps.append(p) for p in self.DB.cursor]
         for ind,f in zip(range(len(Fields)),Fields):#f in Fields:
            #self.tmp.__setattr__(f,np.array([line[f] for line in self.DB.cursor]))
            self.tmp.__setattr__(f,np.array([p[ind] for p in ps]))
    
    def plottrajectory(self,ax,p):
        ##### retrieve #######
        self.part(p)
        self.Part.fields(['Y','X','t'])
        i=self.Part.data.t.argsort()
        ax.plot(self.Part.data.t[i],self.Part.data.Y[i])
        ax.plot(self.Part.data.t[i],self.Part.data.X[i])
    
        
    def detect_islands(self,display=False,displaypic=0,save=True):
        ### Parameters of the simulation ###
        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"
        if os.path.isdir(self.infos.path_rawdata+"DetectionCheck") is False: os.mkdir(self.infos.path_rawdata+"DetectionCheck")
        pathcheck=self.infos.path_rawdata+"DetectionCheck/"+self.infos.suffix+"/"
        if os.path.isdir(pathcheck) is False: os.mkdir(pathcheck)
                
        pics=os.listdir(path+'Frames')        
        if display is True:
            fig=plt.figure(1,figsize=(10,10))
            ax1=fig.add_subplot(111)
            fig2=plt.figure(2,figsize=(10,10))
            ax2=fig2.add_subplot(111)
            #ax2=fig.add_subplot(122)
        else:
            ax1=None;ax2=None        
        num=[]
        numstr=[]     #@@@ CHANGE MADE   26-07-2019 
        for ff in os.listdir(path+'Frames'):
            if ".tiff" in ff:
                numstr+=[ff[-9:-5]]        #@@@ CHANGE MADE   26-07-2019
                num+=[int(ff[-9:-5])]      #@@@ CHANGE MADE   26-07-2019
                #print(num)
        
        picnamestart = os.listdir(path+'Frames')[0][:-9] #@@@ CHANGE MADE   26-07-2019
        num=np.array(num)
        num.sort()
        pics=num.tolist()
        print(pics)
        self.load(dat='timeseries');
        self.timeseries.t=np.zeros(len(pics))
        self.timeseries.pos2=np.zeros(len(pics))
        tt=0
        self.imagedata.islands=[[]]*len(pics)
        if display is True:pics=[pics[displaypic]]
        #print(pics)
        indpic=-1      #@@@ CHANGE MADE   26-07-2019
        for jj in pics:
            indpic=indpic+1     #@@@ CHANGE MADE   26-07-2019
            file_in= picnamestart+numstr[indpic]+'.tiff' #@@@ CHANGE MADE   26-07-2019 from: file_in= "frame_%04d.tiff" % (jj)
            file_out=picnamestart+numstr[indpic]+'.jpg'            #@@@ CHANGE MADE   26-07-2019 from:   file_out= "frame_%04d.jpg" % (jj)
            #Load image
            print(path +'Frames/'+ file_in)
            originalImage = cv2.imread(path + file_in)
            #img = cv2.imread(path+'Frames/' + file_in,0)
            try:
            #if originalImage is not None: #Does the image exist?
                #///////////////////// INITIAL STEP ///////////////////////
                ######## Load ####
                img = cv2.imread(path+'Frames/' + file_in,0)#load image as grayscale image
                #back = cv2.imread(path+'Background/' + file_in,0)
                CLONE = np.copy(img)      # to draw check picture  
                CLONE3 = np.copy(originalImage)
                ######## Rotate and Crop ####
                x0,deltax,y0,deltay=self.infos.BBox
                imgcrop = np.copy(img[y0:y0+deltay, x0:x0+deltax])#Crop image)
                #backcrop = np.copy(back[y0:y0+deltay, x0:x0+deltax])#Crop image)
                del originalImage
                #imgray = imgcrop.astype('int16') 
                #backgray = backcrop.astype('int16') 
                #img2=(imgray-backgray)/2+125
                #img2 = img2.astype('uint8')     
                img2=imgcrop
                if display is True:
                    ax1.imshow(img2)
                (tmp,ThresholdedImage)= cv2.threshold(img2,self.infos.ths,255,cv2.THRESH_BINARY_INV)
                #(tmp,ThresholdedImage)= cv2.threshold(imgcrop,255-self.infos.ths2,255-self.infos.ths1,cv2.THRESH_BINARY_INV)
                del tmp
                if display is True:
                    ax2.imshow(ThresholdedImage)
                #stop
                #############################################
                ################### STEP 1 ##################
                #############################################
                ######### Erode and Dilate : open holes #####
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.erode(ThresholdedImage,kernel)
                kernel = np.ones((5,5),np.uint8)
                ThresholdedImage=cv2.dilate(ThresholdedImage,kernel)
                #ax2.imshow(ThresholdedImage)
                ######## Extract contours ####    
                #ThresholdedImage = ThresholdedImage.astype('uint8') 
                #ax2.imshow(img)                
                im2, contours, hierarchy = cv2.findContours(ThresholdedImage,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

                #Loop over contours to save position and area of islands
                islands=[]
                for cc in contours:
                    mom=cv2.moments(cc)
                    x=mom["m10"]/mom["m00"]
                    y=mom["m01"]/mom["m00"]
                    area=mom["m00"]
                    if area > self.infos.minArea and area<self.infos.maxArea:
                        #if display is True:                        
                        #    ax2.plot(x,y,marker=".",c="w",markersize=12,markeredgewidth=6)
                        cv2.drawContours(imgcrop,[cc],-1,(255,255,255),-1)
                        islands+=[((x,y),area)]                
                self.imagedata.islands[tt]=islands
                cv2.imwrite(pathcheck+file_out,imgcrop)
                del islands,ThresholdedImage,imgcrop,CLONE,CLONE3
                tt+=1
                    ######## Analyze #######
            except: 
                traceback.print_exc()                
                pass
        self.save(dat="imagedata")
    
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
            #box=ellipse2square([[self.Time.data.X[e]*scaling+offset[0],self.Time.data.Y[e]*scaling+offset[1]],[],self.Time.data.phi[e]*np.pi/180],squaresize=np.sqrt(2)/2*scaling,inte=False)
            ax.scatter(self.Time.data.X[e]*scaling+offset[0],self.Time.data.Y[e]*scaling+offset[0])
            
    def track_islands(self):
        self.load(dat="imagedata")
        #Writes a preliminary text file, to be used by the tracking program afterwards
        pathfile=self.path + self.name + '/Files/imagedata_raw' + str(self.sub) + '.txt'
        text_file = open(pathfile, 'wb')
        #Header of the text file
        text_file.writelines([b't\t p\t x0\t y0\t area\n'])
        #start = time.time()
        for jj in range(len(self.imagedata.islands)):
            for kk in range(len(self.imagedata.islands[jj])):
                text_file.write(b"%i\t%i\t%f\t%f\t%f\n" % (jj,kk,self.imagedata.islands[jj][kk][0][0],self.imagedata.islands[jj][kk][0][1],self.imagedata.islands[jj][kk][1]))
        text_file.close()
        #Opens the text file and do the tracking
        self.imagedata.tpos=np.loadtxt(pathfile,skiprows=1)
        self.imagedata.P=tracking(self.imagedata.tpos,usesparse=True,display=True,mem=2,delta=self.infos.tracking_maxdisp)

        self.save(dat="imagedata")
        
    def correct_displacement(self):
        self.load(dat='timeseries')
        x=np.linspace(0,2048,200)
        f1=interp1d(np.polyval(self.infos.distortion,x),x)
        self.timeseries.pos2_corrected=f1(self.infos.BBox[0]+self.timeseries.pos2)*self.infos.pix2mm
        self.save(dat='timeseries')





    def plot_ellipses_orderparam(self,ax1,jj,t):
        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"
        picnamestart = os.listdir(path+'Frames')[0][:-9] #@@@ CHANGE MADE   30-07-2019
        file_in= picnamestart+"%04d.tiff" % (jj)  #@@@ CHANGE MADE   30-07-2019: file_in= "frame_%04d.tiff" % (jj)
        img = cv2.imread(path+'Frames/' + file_in,0)
        ax1.imshow(img,cmap="gray")                        
        scaling=1/self.infos.pix2mm*self.infos.square
        sizes=0.8
        offset=[self.infos.BBox[0],self.infos.BBox[2]]
        self.time(t)
        self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        
        req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (t=%s) ORDER BY e;' % (self.sub,self.Time.t)
        self.fetch(req,['e','X','Y','a','b','angle','t'])            
        t=np.copy(self.tmp.t);elab=np.copy(self.tmp.e);
        X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);a=np.copy(self.tmp.a);b=np.copy(self.tmp.b);angle=np.copy(self.tmp.angle)
        
        Omega=(1-a/b)*np.cos(2*(angle-np.pi/4))
        for e in range(len(X)):
            ix=np.where(eref==elab[e])
            el=Ellipse(xy=[X[e]*scaling+offset[0],Y[e]*scaling+offset[1]],
                       width=a[e]*scaling*sizes,height=b[e]*scaling*sizes,angle=angle[e]*180/np.pi)
            ax1.add_artist(el)
            #el.set_alpha(0.5)
            ialpha=np.abs(Omega[e])*5
            if ialpha >=1: ialpha=1.0
            
            #ialpha=1-ialpha
            #ialpha=1
            #ialpha=1.
            el.set_alpha(ialpha)
            #icol=(Omega[e]+1)/2#icol=(Omega[e]*(-1)**(nxref[ix][0]+nyref[ix][0])+1)/2.
            icol=(Omega[e]*(-1)**(nyref[ix][0])+1)/2.
            el.set_facecolor(plt.cm.bwr(icol))#plt.cm.PuOr(icol))                
            #print('Om ='+str(Omega[e])+' , icol ='+str(icol))
            el.set_edgecolor(None)
            ax1.add_artist(el)

    












        
    def plot_ellipses_orderparam_BACKUP(self,ax1,jj,t):
        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"
        picnamestart = os.listdir(path+'Frames')[0][:-9] #@@@ CHANGE MADE   30-07-2019
        file_in= picnamestart+"%04d.tiff" % (jj)  #@@@ CHANGE MADE   30-07-2019: file_in= "frame_%04d.tiff" % (jj)
        img = cv2.imread(path+'Frames/' + file_in,0)
        ax1.imshow(img,cmap="gray")                        
        scaling=1/self.infos.pix2mm*self.infos.square
        sizes=0.8
        offset=[self.infos.BBox[0],self.infos.BBox[2]]
        self.time(t)
        self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (t=%s) ORDER BY e;' % (self.sub,self.Time.t)
        self.fetch(req,['e','X','Y','a','b','angle','t'])            
        t=np.copy(self.tmp.t);elab=np.copy(self.tmp.e);
        X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);a=np.copy(self.tmp.a);b=np.copy(self.tmp.b);angle=np.copy(self.tmp.angle)
        
        Omega=(1-a/b)*np.cos(2*angle)
        for e in range(len(X)):
            ix=np.where(eref==elab[e])
            el=Ellipse(xy=[X[e]*scaling+offset[0],Y[e]*scaling+offset[1]],
                       width=a[e]*scaling*sizes,height=b[e]*scaling*sizes,angle=angle[e])
            ax1.add_artist(el)
            el.set_alpha(0.5)
            icol=(Omega[e]*(-1)**(nxref[ix][0]+nyref[ix][0])+1)/2.
            el.set_facecolor(plt.cm.bwr(icol))#plt.cm.PuOr(icol))
            el.set_edgecolor(None)


    def plot_orderparam_symaxes(self,ax1,jj,t):



        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"

        self.time(t)
        self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        #req='SELECT s,nx,ny FROM SquareProp%s ORDER BY s;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (t=%s) ORDER BY e;' % (self.sub,self.Time.t)
        #req='SELECT s,X,Y,a,b,angle,t FROM Square%s WHERE (t=%s) ORDER BY s;' % (self.sub,self.Time.t)
        self.fetch(req,['e','X','Y','a','b','angle','t'])            
        t=np.copy(self.tmp.t);elab=np.copy(self.tmp.e);
        X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);a=np.copy(self.tmp.a);b=np.copy(self.tmp.b);angle=np.copy(self.tmp.angle)
        
        
        
        #Omega=(1-a/b)*np.cos(2*angle)
        Omega=(1-a/b)*np.cos(2*(angle-np.pi/4))

        nxcent = np.int(np.average(np.unique(nxref)))
        nycent = np.int(np.average(np.unique(nyref)))
        
        #nyref_inv = np.array(len(X)*[int(np.max(nyref)+1)])-nyref
        
        Omxcent = []
        Omycent = []
        
        Omxcent_args = []
        Omycent_args = []
        output_dic={}

        Omxcent = Omega[np.argwhere(nxref==nxcent)]
        Omycent = Omega[np.argwhere(nyref==nycent)]


        
        
        for ind in range(len(Omxcent)):
            Omxcent[ind]=Omxcent[ind][0]
        for ind in range(len(Omycent)):
            Omycent[ind]=Omycent[ind][0]


        
        ax1.plot(Omxcent,label='y-axis')
        ax1.plot(Omycent,label='x-axis')
        
        ax1.set_xlabel('position')
        ax1.set_ylabel('polarisation P')
        if  np.max(np.abs(np.array(Omxcent+Omycent)))<0.1:
            ax1.set_ylim(-0.1,0.1)                
            ax1.set_yticks([-0.1,0,0.1])
        elif  np.max(np.abs(np.array(Omxcent+Omycent)))> 0.1 and np.max(np.abs(np.array(Omxcent+Omycent)))< 0.5 :
            ax1.set_ylim(-0.5,0.5)                
            ax1.set_yticks([-0.5,0,0.5])
        elif  np.max(np.abs(np.array(Omxcent+Omycent)))> 0.5:
            ax1.set_ylim(-1,1)                
            ax1.set_yticks([-1,0,1])        
        ax1.xaxis.set_label_coords(0.5, -0.1)
        ax1.yaxis.set_label_coords(-0.13, 0.5)
        
        
        Nx = len(Omxcent)
        Ny = len(Omycent)
        N = np.max([Nx,Ny])
        emp_ar = np.full(N,np.nan)
        
        pos_x = np.array(range(Nx),dtype=float)
        pos_y = np.array(range(Ny),dtype=float)

        
        output_dic['Position x']=np.concatenate((    pos_x   ,emp_ar))[:N]
        output_dic['Omega x, centre']=np.concatenate((    Omxcent[:,0]   ,emp_ar))[:N]
        output_dic['Position y']=np.concatenate((    pos_y   ,emp_ar))[:N]
        output_dic['Omega y, centre']=np.concatenate((    Omycent[:,0]   ,emp_ar))[:N]

        output_DF = pd.DataFrame(output_dic)
            
        return output_DF
 


        


    def plot_orderparam_xyavg(self,ax1,jj,t):



        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"

        self.time(t)
        self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        #req='SELECT s,nx,ny FROM SquareProp%s ORDER BY s;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (t=%s) ORDER BY e;' % (self.sub,self.Time.t)
        #req='SELECT s,X,Y,a,b,angle,t FROM Square%s WHERE (t=%s) ORDER BY s;' % (self.sub,self.Time.t)
        self.fetch(req,['e','X','Y','a','b','angle','t'])            
        t=np.copy(self.tmp.t);elab=np.copy(self.tmp.e);
        X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);a=np.copy(self.tmp.a);b=np.copy(self.tmp.b);angle=np.copy(self.tmp.angle)
        
        
        
        #Omega=(1-a/b)*np.cos(2*angle)
        Omega_reg=(1-a/b)*np.cos(2*(angle-np.pi/4))
        
        eveny = np.argwhere(nyref/2==np.floor(nyref/2))
        oddy = np.argwhere(nyref/2!=np.floor(nyref/2))
        
        Omega = 1.*Omega_reg
        Omega[eveny]=-1.*Omega[eveny]
        nxunique = np.unique(nxref)
        nyunique = np.unique(nyref)
        
        Omx = np.zeros(len(nxunique))
        Omy = np.zeros(len(nyunique))
        
        
        for ind in range(len(Omx)):
            Omx[ind] = np.average(Omega[np.argwhere(nxref==nxunique[ind])])
        for ind in range(len(Omy)):
            Omy[ind] = np.average(Omega[np.argwhere(nyref==nyunique[ind])])
            
        Omx =list(Omx)
        Omy =list(Omy)
        
        ax1.plot(Omx,label='x-direction')
        ax1.plot(Omy,label='y-direction')
        
        
        ax1.set_xlabel('position')
        ax1.set_ylabel('polarisation P')
        if  np.max(np.abs(np.array(Omx+Omy)))<0.1:
            ax1.set_ylim(-0.1,0.1)                
            ax1.set_yticks([-0.1,0,0.1])
        elif  np.max(np.abs(np.array(Omx+Omy)))> 0.1 and np.max(np.abs(np.array(Omx+Omy)))< 0.5 :
            ax1.set_ylim(-0.5,0.5)                
            ax1.set_yticks([-0.5,0,0.5])
        elif  np.max(np.abs(np.array(Omx+Omy)))> 0.5:
            ax1.set_ylim(-1,1)                
            ax1.set_yticks([-1,0,1])        
        ax1.xaxis.set_label_coords(0.5, -0.1)
        ax1.yaxis.set_label_coords(-0.1, 0.5)
        
        
        Nx = len(Omx)
        Ny = len(Omy)
        N = np.max([Nx,Ny])
        
        pos_x = np.array(range(Nx),dtype=float)
        pos_y = np.array(range(Ny),dtype=float)
        
        emp_ar = np.full(N,np.nan)
        output_dic={}
        output_dic['Position x']=np.concatenate((    pos_x   ,emp_ar))[:N]
        output_dic['Omega x, centre']=np.concatenate((    np.array(Omx),emp_ar))[:N]
        output_dic['Position y']=np.concatenate((    pos_y,emp_ar))[:N]
        output_dic['Omega y, centre']=np.concatenate((    np.array(Omy),emp_ar))[:N]
        
        
        output_DF = pd.DataFrame(output_dic)
        
        return output_DF


        
    def plot_orderparam_diag(self,ax1,jj,t):



        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"

        self.time(t)
        self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (t=%s) ORDER BY e;' % (self.sub,self.Time.t)
        self.fetch(req,['e','X','Y','a','b','angle','t'])            
        t=np.copy(self.tmp.t);elab=np.copy(self.tmp.e);
        X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);a=np.copy(self.tmp.a);b=np.copy(self.tmp.b);angle=np.copy(self.tmp.angle)
        
        Omega=(1-a/b)*np.cos(2*angle)
        
        nyref_inv = np.array(len(X)*[int(np.max(nyref)+1)])-nyref
        
        Omdiag1 = []
        Omdiag2 = []
        
        
        for e in range(len(X)):
            if nxref[e]==nyref[e]:
                Omdiag1+=[Omega[e]]
            elif nxref[e]==nyref_inv[e]:
                Omdiag2+=[Omega[e]]
        
        ax1.plot(Omdiag1)
        ax1.plot(Omdiag2)
        ax1.set_xlabel('x position')
        ax1.set_ylabel('polarisation P')
        if  np.max(np.abs(np.array(Omdiag1+Omdiag2)))<0.1:
            ax1.set_ylim(-0.1,0.1)                
            ax1.set_yticks([-0.1,0,0.1])
        elif  np.max(np.abs(np.array(Omdiag1+Omdiag2)))> 0.1 and np.max(np.abs(np.array(Omdiag1+Omdiag2)))< 0.5 :
            ax1.set_ylim(-0.5,0.5)                
            ax1.set_yticks([-0.5,0,0.5])
        elif  np.max(np.abs(np.array(Omdiag1+Omdiag2)))> 0.5:
            ax1.set_ylim(-1,1)                
            ax1.set_yticks([-1,0,1])        
        ax1.xaxis.set_label_coords(0.5, -0.1)
        ax1.yaxis.set_label_coords(-0.13, 0.5)

        
        """
        for e in range(len(X)):
            ix=np.where(eref==elab[e])
            el=Ellipse(xy=[X[e]*scaling+offset[0],Y[e]*scaling+offset[1]],
                       width=a[e]*scaling*sizes,height=b[e]*scaling*sizes,angle=angle[e])
            ax1.add_artist(el)
            el.set_alpha(0.5)
            icol=(Omega[e]*(-1)**(nxref[ix][0]+nyref[ix][0])+1)/2.
            el.set_facecolor(plt.cm.bwr(icol))#plt.cm.PuOr(icol))
            el.set_edgecolor(None)
        """






    def plot_Poisson_el(self,ax1,ax2,nxy):



        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"
        
        #self.time(t)
        #self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        #req='SELECT s,nx,ny FROM SquareProp%s ORDER BY s;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        xycomb=[    [nxy[0],nxy[0+2]],     [nxy[0],nxy[1+2]],     [nxy[1],nxy[0+2]],     [nxy[1],nxy[1+2]]    ]
        
        evec = np.array(4*[0],dtype=int)
        xvec = np.array(4*[0],dtype=int)
        yvec = np.array(4*[0],dtype=int)
        for ind in range(len(nxy)):
            indx = np.where(nxref==xycomb[ind][0] )
            indy = np.where(nyref==xycomb[ind][1])
            indxy=np.intersect1d(indx,indy)
            evec[ind]=eref[indxy[0]]
            xvec[ind]=nxref[indxy[0]]
            yvec[ind]=nyref[indxy[0]]
            
            
        t=4*[None]
        elab=4*[None]
        X=4*[None]
        Y=4*[None]
        a=4*[None]
        b=4*[None]
        angle=4*[None]
        t=4*[None]
        for ind in range(len(nxy)):
            
            e=evec[ind]
            req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (e=%e) ORDER BY t;'  % (self.sub,e)#,self.Time.t)
            
            
            #req='SELECT s,X,Y,a,b,angle,t FROM Square%s WHERE (t=%s) ORDER BY s;' % (self.sub,self.Time.t)
            self.fetch(req,['e','X','Y','a','b','angle','t'])    
            
            t[ind]=np.copy(self.tmp.t);elab[ind]=np.copy(self.tmp.e);
            X[ind]=np.copy(self.tmp.X);Y[ind]=np.copy(self.tmp.Y);a[ind]=np.copy(self.tmp.a);b[ind]=np.copy(self.tmp.b);angle[ind]=np.copy(self.tmp.angle)

        lenan = []
        for ind in range(len(nxy)):
            lenan+=[len(t[ind])]
        lenan=np.min(np.array(lenan,dtype=int))
        
        for ind in range(len(nxy)):
            t[ind]=t[ind][:lenan]
            elab[ind]=elab[ind][:lenan]
            X[ind]=X[ind][:lenan]
            Y[ind]=Y[ind][:lenan]
            a[ind]=a[ind][:lenan]
            b[ind]=b[ind][:lenan]
            angle[ind]=angle[ind][:lenan]
            
        
        LY = (Y[1]+Y[3])/2-(Y[0]+Y[2])/2
        LX = (X[2]+X[3])/2-(X[0]+X[1])/2
        
        epsY = (LY-LY[0])/LY[0] 
        epsX = (LX-LX[0])/LX[0] 
        nu = -epsX/epsY
        
        
        comp_start= np.where(epsY<0.05*np.min(epsY))[0][0]        
        comp_end=np.where(epsY<0.95*np.min(epsY))[0][0]
        
        
        epsY_ad = np.concatenate((np.array([0.]),epsY[comp_start:comp_end]))
        epsX_ad = np.concatenate((np.array([0.]),epsX[comp_start:comp_end]))
        #nu_ad = -epsX_ad/epsY_ad
        
        steps_disp = 20
        locs = np.linspace(0,len(epsY_ad)-1,steps_disp,dtype=int)
        epsY_ad=epsY_ad[locs]
        epsX_ad=epsX_ad[locs]
        nu_ad = -epsX_ad/epsY_ad
        
        
        
        depsY = epsY[1:]-epsY[:-1]
        depsX = epsX[1:]-epsX[:-1]
        dnu=-depsX/depsY
        
        
        
        
        depsY_ad = epsY_ad[1:]-epsY_ad[:-1]
        depsX_ad = epsX_ad[1:]-epsX_ad[:-1]
        dnu_ad=-depsX_ad/depsY_ad
        
        ax1.plot(-epsY_ad[1:],nu_ad[1:])
        ax2.plot(-epsY_ad[1:],dnu_ad)
        
        
        ax1.set_xlabel('$-\epsilon_y$')
        ax1.set_ylabel('$ \\nu $')
        ax1.xaxis.set_label_coords(0.5, -0.1)
        ax1.yaxis.set_label_coords(-0.15, 0.5)
        
        ax2.set_xlabel('$-\epsilon_y$')
        ax2.set_ylabel('$ \\nu $')
        ax2.xaxis.set_label_coords(0.5, -0.1)
        ax2.yaxis.set_label_coords(-0.15, 0.5)
                

        
        
        output_dic={}
        output_dic['- epsilon_y']=-epsY_ad[1:]
        output_dic['nu AveragePoisson']=nu_ad[1:]
        output_dic['nu RelativePoisson']=dnu_ad
        output_DF = pd.DataFrame(output_dic)
        
        return output_DF
        
        
        
        
        
        """
        ax1.set_xlabel('position')
        ax1.set_ylabel('polarisation P')
        if  np.max(np.abs(np.array(Omx+Omy)))<0.1:
            ax1.set_ylim(-0.1,0.1)                
            ax1.set_yticks([-0.1,0,0.1])
        elif  np.max(np.abs(np.array(Omx+Omy)))> 0.1 and np.max(np.abs(np.array(Omx+Omy)))< 0.5 :
            ax1.set_ylim(-0.5,0.5)                
            ax1.set_yticks([-0.5,0,0.5])
        elif  np.max(np.abs(np.array(Omx+Omy)))> 0.5:
            ax1.set_ylim(-1,1)                
            ax1.set_yticks([-1,0,1])        
        ax1.xaxis.set_label_coords(0.5, -0.1)
        ax1.yaxis.set_label_coords(-0.13, 0.5)

        """
        
        

    def plot_Poisson_el_avg(self,ax1,ax2,ax3,nxy,epsYlim):



        path=self.infos.path_rawdata+"Rawdata/"+self.infos.suffix+"/"
        
        #self.time(t)
        #self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s ORDER BY e;' % (self.sub)
        #req='SELECT s,nx,ny FROM SquareProp%s ORDER BY s;' % (self.sub)
        self.fetch(req,['e','nx','ny'])
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        #xycomb=[    [nxy[0],nxy[0+2]],     [nxy[0],nxy[1+2]],     [nxy[1],nxy[0+2]],     [nxy[1],nxy[1+2]]    ]
        
        
        
        
        evec = 4*[0]
        xvec = 4*[0]
        yvec = 4*[0]
        
        for ind in [0,1]:
            evec[ind] = eref[np.where(nxref==nxy[ind])[0]]
            xvec[ind] = nxref[np.where(nxref==nxy[ind])[0]]
            yvec[ind] = nyref[np.where(nxref==nxy[ind])[0]]            
        for ind in [2,3]:
            evec[ind] = eref[np.where(nyref==nxy[ind])[0]]
            xvec[ind] = nxref[np.where(nyref==nxy[ind])[0]]
            yvec[ind] = nyref[np.where(nyref==nxy[ind])[0]]
            

            
            
            
        """
        t=4*[{}]
        elab=4*[{}]
        X=4*[{}]
        Y=4*[{}]
        a=4*[{}]
        b=4*[{}]
        angle=4*[{}]
        
        """
        t=4*[None]
        elab=4*[None]
        X=4*[None]
        Y=4*[None]
        a=4*[None]
        b=4*[None]
        angle=4*[None]
        
        for var in [t,elab,X,Y,a,b,angle]:
            for ind in range(4):
                var[ind]=len(evec[ind])*[None]
        
        
        for ind in range(len(nxy)):
            for ind2 in range(len(evec[ind])):
            
                e=copy.deepcopy(evec[ind][ind2])
                
                req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (e=%e) ORDER BY t;'  % (self.sub,e);                self.fetch(req,['e','X','Y','a','b','angle','t'])    
                
                
                t[ind][ind2]=copy.deepcopy(self.tmp.t);           elab[ind][ind2]=copy.deepcopy(self.tmp.e);
                X[ind][ind2]=copy.deepcopy(self.tmp.X);       Y[ind][ind2]=copy.deepcopy(self.tmp.Y);       
                a[ind][ind2]=copy.deepcopy(self.tmp.a);       b[ind][ind2]=copy.deepcopy(self.tmp.b);   angle[ind][ind2]=copy.deepcopy(self.tmp.angle)
                """
                t[ind][ind2]=np.copy(self.tmp.t);           elab[ind][ind2]=np.copy(self.tmp.e);
                X[ind][ind2]=np.copy(self.tmp.X);       Y[ind][ind2]=np.copy(self.tmp.Y);       
                a[ind][ind2]=np.copy(self.tmp.a);       b[ind][ind2]=np.copy(self.tmp.b);   angle[ind][ind2]=np.copy(self.tmp.angle)
                """
                #print(X[ind][ind2][:10])


        lenan = []
        for ind in range(len(nxy)):
            for ind2 in range(len(evec[ind])):
                
                lenan+=[len(t[ind][ind2])]
        lenan=np.min(np.array(lenan,dtype=int))
        
    
        for ind in range(len(nxy)):
            for ind2 in range(len(evec[ind])):
                t[ind][ind2]        =t[ind][ind2][:lenan]
                elab[ind][ind2]     =elab[ind][ind2][:lenan]
                X[ind][ind2]        =X[ind][ind2][:lenan]
                Y[ind][ind2]        =Y[ind][ind2][:lenan]
                a[ind][ind2]        =a[ind][ind2][:lenan]
                b[ind][ind2]        =b[ind][ind2][:lenan]
                angle[ind][ind2]    =angle[ind][ind2][:lenan]
        
        
        Xavg=4*[None]
        Yavg=4*[None]
        
        for ind in range(4):
            Xavg[ind]=np.zeros(lenan)
            Yavg[ind]=np.zeros(lenan)
            for ind2 in range(len(X[ind])):
                Xavg[ind]+=X[ind][ind2]/len(X[ind])
                Yavg[ind]+=Y[ind][ind2]/len(Y[ind])
                
        
        LY = Yavg[3]-Yavg[2]
        LX = Xavg[1]-Xavg[0]
        
        epsY = (LY-LY[0])/LY[0] 
        epsX = (LX-LX[0])/LX[0] 
        
        
        
        interpf = interp1d(epsY,epsX,'linear')
        
        ninterp = 1000
        n_SG = 4
        window_SG = 2*int(ninterp/n_SG/2)+1
        order_SG=3
        
        
        epsY_interp = np.linspace(0.,epsYlim,ninterp)
        epsX_interp= interpf(epsY_interp)

        epsX_SG = savgol_filter(epsX_interp,window_SG,order_SG)
        
        
        
        
        
        
        DepsX = -1.*epsX_SG[0]
        epsX_SG+=DepsX
        epsX_interp+=DepsX
        epsX+=DepsX
        
        n_deps = 10
        
        epsY_n = np.linspace(0.,min(epsY_interp),n_deps+1)
        
        
        interpfn = interp1d(epsY_interp,epsX_SG,'linear')
        epsX_n = interpfn(epsY_n)
        
        
        
        
        epsY_final=copy.deepcopy(epsY_interp)
        epsX_final=copy.deepcopy(epsX_SG)
        
        
        

        ax3.plot(-epsY,epsX,label='Data')
        ax3.plot(-epsY_interp,epsX_SG, label='Savitzky-Golay')
        ax3.set_xlabel('-$\epsilon_{y}$')
        ax3.set_ylabel('$\epsilon_{x}$')
        ax3.set_xlim(0.,-epsYlim)
        
        
        
        nu = -epsX_final[1:]/epsY_final[1:]
        epsY_nu = epsY_final[1:]
        nu_SG = savgol_filter(nu,window_SG,order_SG)

        
        
        
        epsY_dnu = epsY_n[1:]
        
        depsY = epsY_n[1:]-epsY_n[:-1]
        depsX = epsX_n[1:]-epsX_n[:-1]
        dnu=-depsX/depsY
        #dnu_SG = savgol_filter(dnu,window_SG,order_SG)

        
        
        ax1.plot(-epsY_nu,nu)
        ax1.plot(-epsY_nu,nu_SG)
        ax1.set_xlabel('$-\epsilon_y$')
        ax1.set_ylabel('$ \\nu $')
        ax1.xaxis.set_label_coords(0.5, -0.1)
        ax1.yaxis.set_label_coords(-0.15, 0.5)
        ax1.set_xlim(0.,-epsYlim)


        ax2.plot(-epsY_dnu,dnu,'--',marker='o',markersize=20, markerfacecolor='none')        
        ax2.set_xlabel('$-\epsilon_y$')
        ax2.set_ylabel('$ \\nu $')
        ax2.xaxis.set_label_coords(0.5, -0.1)
        ax2.yaxis.set_label_coords(-0.15, 0.5)
        ax2.set_xlim(0.,0.06)#-epsYlim)
        ax2.set_ylim(-0.6,0.6)

        
        output_dic={}
        output_dic['- epsilon_y RelativePoisson_avgside']=-epsY_dnu
        output_dic['nu RelativePoisson_avgside']=dnu        
        output_DF = pd.DataFrame(output_dic)
        
        return output_DF
        
        














class CreateDB(Focus):
    
    def __init__(self,name,sub):
        self.path=path_sys[_platform]+'Experiments/'
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
                t=pp[:,0];X=pp[:,1]*self.infos.pix2mm/self.infos.square;Y=pp[:,2]*self.infos.pix2mm/self.infos.square;
                area=pp[:,3]*self.infos.pix2mm**2/self.infos.square**2;
                tt=0
                #print(len(t))
                for i in range(len(t)):
                    fid.write('%d;%d;%f;%f;%f\n'% (p,t[i],X[i],Y[i],area[i]))
                    fid2.write('%d\n'% (p))
                    tt+=1
            except:
                print(p)
                t=pp[0];X=pp[1]*self.infos.pix2mm/self.infos.square;Y=pp[2]*self.infos.pix2mm/self.infos.square;
                area=pp[3]*self.infos.pix2mm**2/self.infos.square**2;
                fid.write('%d;%d;%f;%f;%f\n'% (p,t,X,Y,area))
                fid2.write('%d\n'% (p))
        fid.close()
        fid2.close()
       
        
        
        
        
        
        
        
        
        
    def writetxtSQUARE(self,ax1=None):          
        #retrieve data
        req='SELECT DISTINCT nx FROM Prop%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx']) 
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM Prop%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
        #initialise folders
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Square.txt','w')
        fid2=open(tmpdir + '/SquareProp.txt','w')
        
        #loop over columns
        s=0
        for nx in nxs:
            for ny in nys:
                req='SELECT p FROM Prop%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
                #print(req)
                self.fetch(req,['p'])
                ps=self.tmp.p
                square=False
            
                #try:
                """
                if len(ps) == 4: 
                    #calculate_cm_orientation(self.tmp.p)
                    p=ps[0]            
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    
                    Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    
                    #Xvec=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                    #Yvec=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                    Xvec=[list(np.copy(self.Part.data.X)[:self.infos.tmax])]
                    Yvec=[list(np.copy(self.Part.data.Y)[:self.infos.tmax])]
                    
                    Phicm=np.zeros(self.infos.tmax)
                    for p in ps[1:]:
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                        Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                        Xvec+=[list(np.copy(self.Part.data.X)[:self.infos.tmax])]
                        Yvec+=[list(np.copy(self.Part.data.Y)[:self.infos.tmax])]
                    
                    Xcm=Xcm/4
                    Ycm=Ycm/4
                    Xvec=np.array(Xvec)
                    Yvec=np.array(Yvec)
                    
                    
                    Phicm=np.zeros(len(Xcm))
                    for ind in range(len(Phicm)):
                        
                        Phicm[ind] = np.arctan(      (np.average(Yvec[np.argwhere(Xvec[:,ind]>Xcm[ind]),ind]) -Ycm[ind])   /(np.average(Xvec[np.argwhere(Xvec[:,ind]>Xcm[ind]),ind]) -Xcm[ind]))
                    #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                    #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                    square=True
                                    
                """




                    
                if len(ps) > 2.9: 
                    #calculate_cm_orientation(self.tmp.p)
                    tcm=[]
                    Xcm=[]
                    Ycm=[]
                    ltrack=[]#np.zeros(len(ps),dtype=int)
                    for ind in range(len(ps)):
                        p=ps[ind]            
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        tn = np.copy(self.Part.data.t)[:self.infos.tmax];

                        if tn[0]==0:    
                            tcm+=[np.copy(self.Part.data.t)[:self.infos.tmax]]
                            Xcm+=[np.copy(self.Part.data.X)[:self.infos.tmax]]
                            Ycm+=[np.copy(self.Part.data.Y)[:self.infos.tmax]]
                            ltrack+=[len(tcm[ind])]
                    if len(tcm)>2.9:
                        ltrack=np.array(ltrack,dtype=int)
                        tmaxtrack=np.sort(ltrack)[-3]
                        ps_inds=np.argsort(ltrack)[-3:]
                        
                        indcomb1 = [ps_inds[1],ps_inds[2],ps_inds[2]]
                        indcomb2 = [ps_inds[0],ps_inds[0],ps_inds[1]]
                        dvec0 = np.zeros(3)
                        for ind in range(3):
                            ind1 = indcomb1[ind]
                            ind2 = indcomb2[ind]
                            
                            dvec0[ind] = np.sqrt(    (Xcm[ind1][0]-Xcm[ind2][0])**2   +  (Ycm[ind1][0]-Ycm[ind2][0])**2)
                            
                        combmax=np.argmax(dvec0)
                        indsdiag = [indcomb1[combmax],indcomb2[combmax]]
                        
                        Xmat = np.zeros((2,tmaxtrack))
                        Ymat = np.zeros((2,tmaxtrack))
                        for ind in range(2):
                            Xmat[ind,:] = 1.*Xcm[indsdiag[ind]][:tmaxtrack]
                            Ymat[ind,:] = 1.*Ycm[indsdiag[ind]][:tmaxtrack]
                        
                        Xcm = (Xmat[0,:]+Xmat[1,:])/2
                        Ycm = (Ymat[0,:]+Ymat[1,:])/2
                        
                        dX=Xmat[1,:]-Xmat[0,:]
                        dY=Ymat[1,:]-Ymat[0,:]
                            
                        
                        Phicm = np.arctan(dY/dX)-np.arctan(dY[0]/dX[0])      
                        tcm = np.array(range(tmaxtrack),dtype=int)
    
                    
                        #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                        #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                        square=True
                    
                    
                """                    
                if len(ps) == 3 or len(ps) ==4: 
                    #calculate_cm_orientation(self.tmp.p)
                    tcm=len(ps)*[None]
                    Xcm=len(ps)*[None]
                    Ycm=len(ps)*[None]
                    ltrack=np.zeros(len(ps),dtype=int)
                    for ind in range(len(ps)):
                        p=ps[ind]            
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        tcm[ind]=np.copy(self.Part.data.t)[:self.infos.tmax];
                        Xcm[ind]=np.copy(self.Part.data.X)[:self.infos.tmax];
                        Ycm[ind]=np.copy(self.Part.data.Y)[:self.infos.tmax];
                        ltrack[ind]=len(tcm[ind])
                    tmaxtrack=np.sort(ltrack)[-3]
                    ps_inds=np.argsort(ltrack)[-3:]
                    
                    indcomb1 = [ps_inds[1],ps_inds[2],ps_inds[2]]
                    indcomb2 = [ps_inds[0],ps_inds[0],ps_inds[1]]
                    dvec0 = np.zeros(3)
                    for ind in range(3):
                        ind1 = indcomb1[ind]
                        ind2 = indcomb2[ind]
                        
                        dvec0[ind] = np.sqrt(    (Xcm[ind1][0]-Xcm[ind2][0])**2   +  (Ycm[ind1][0]-Ycm[ind2][0])**2)
                        
                    combmax=np.argmax(dvec0)
                    indsdiag = [indcomb1[combmax],indcomb2[combmax]]
                    
                    Xmat = np.zeros((2,tmaxtrack))
                    Ymat = np.zeros((2,tmaxtrack))
                    for ind in range(2):
                        Xmat[ind,:] = 1.*Xcm[indsdiag[ind]][:tmaxtrack]
                        Ymat[ind,:] = 1.*Ycm[indsdiag[ind]][:tmaxtrack]
                    
                    
                    Xcm = (Xmat[0,:]+Xmat[1,:])/2
                    Ycm = (Ymat[0,:]+Ymat[1,:])/2
                    
                    dX=Xmat[1,:]-Xmat[0,:]
                    dY=Ymat[1,:]-Ymat[0,:]
                        
                    
                    Phicm = np.arctan(dY/dX)-np.arctan(dY[0]/dX[0])      
                    tcm = np.array(range(tmaxtrack),dtype=int)

                    #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                    #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                    square=True
                """

                """
                if len(ps) == 2: 
                    #calculate_cm_orientation(self.tmp.p)
                    p=ps[0]            
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    
                    Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    
                    Xvec=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                    Yvec=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                    
                    Phicm=np.zeros(self.infos.tmax)
                    for p in ps[1:]:
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                        Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                        Xvec+=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                        Yvec+=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                    
                    Xvec=np.array(Xvec)
                    Yvec=np.array(Yvec)
                    
                    
                    dX= Xvec[1]-Xvec[0]
                    dY= Yvec[1]-Yvec[0]
                    
                    if nx ==nxs[0] and ny >nys[0]+0.1 and ny < nys[-1]+0.1:     #left edge
                        edgetype = 0
                        if dY<0:
                            dX=-dX
                            dY=-dY
                        Xvec2 = Xvec-dY
                        Yvec2 = Yvec+dX
                    elif nx ==nxs[-1] and ny >nys[0]+0.1 and ny < nys[-1]+0.1:        #right edge
                        edgetype = 1
                        if dY<0:
                            dX=-dX
                            dY=-dY
                        Xvec2 = Xvec+dY
                        Yvec2 = Yvec-dX                        
                    elif ny ==nys[0] and nx >nxs[0]+0.1 and nx < nxs[-1]+0.1:     #bottom edge
                        edgetype = 2
                        if dX<0:
                            dX=-dX
                            dY=-dY
                        Xvec2 = Xvec+dY
                        Yvec2 = Yvec-dX
                    elif ny ==nys[-1] and nx >nxs[0]+0.1 and nx < nxs[-1]+0.1:        #top edge
                        edgetype = 3
                        if dX<0:
                            dX=-dX
                            dY=-dY
                        Xvec2 = Xvec-dY
                        Yvec2 = Yvec+dX
                    else:
                        print('Error: can only track 2 holes but is not on edge')
                        ERRORlksemnrgv
                        
                    Xvec=np.concatenate((Xvec,Xvec2),axis=0)
                    Yvec=np.concatenate((Yvec,Yvec2),axis=0)
                        
                    Xcm = np.array([np.average(Xvec)])
                    Ycm = np.array([np.average(Yvec)])
                    
                    
                    Phicm = np.arctan(      (np.average(Yvec[np.argwhere(Xvec>Xcm)]) -Ycm)   /(np.average(Xvec[np.argwhere(Xvec>Xcm)]) -Xcm))
                    #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                    #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                    square=True
                """
                        
                        
                    
                """
                if len(ps) == 2: 
                    
                    
                    
                    p=ps[0]            
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    Phicm=np.zeros(self.infos.tmax)
                    p=ps[1]
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    
                    vecX=Xcm-np.copy(self.Part.data.X)[:self.infos.tmax]
                    vecY=Ycm-np.copy(self.Part.data.Y)[:self.infos.tmax]
                    
                    Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                    Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                    
                    offx=-np.abs(vecX)/2
                    offy=np.abs(vecY)/2
                    
                    Xcm=Xcm/2+offx
                    Ycm=Ycm/2+offy
                    
                    Phicm=np.arctan2(vecY,vecX)
                    Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]            
                    square=True
                """
                """
                if nx==4 and ny==3:
                    
                    plt.figure()
                    plt.plot(Xmat[0,:])
                    plt.plot(Xmat[1,:])
                    plt.plot(Ymat[0,:])
                    plt.plot(Ymat[1,:])
                    qwjshfnkmvhbj
                """
                
                #check
                if ax1 is not None and square: 
                    ax1.scatter(Xcm,Ycm,color="r")
                 #Save for database
                if square:
                    fid2.write('%d;%d;%d\n'% (s,nx,ny))
                    for ix in range(len(tcm)):
                        fid.write('%d;%d;%f;%f;%f\n'% (s,tcm[ix],Xcm[ix],Ycm[ix],Phicm[ix]))
                    s+=1
                """
                except:
                    pass
                """

        fid.close()
        fid2.close()
        
        req='DROP TABLE Square' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE SquareProp' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass        
        
        req='CREATE TABLE Square' + str(self.sub) + ' ( s int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL , Phi float NOT NULL , PRIMARY KEY ( t ,  s )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE SquareProp' + str(self.sub) + ' ( s int( 8  )  NOT  NULL ,  nx int( 8  )  NOT  NULL , ny int( 8  )  NOT  NULL , PRIMARY KEY ( s , nx, ny )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Square.txt' INTO TABLE Square" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/SquareProp.txt' INTO TABLE SquareProp" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        
        
        
        
        
        
        
        
        
    def writetxtSQUARE_Backup_v2(self,ax1=None):          
        #retrieve data
        req='SELECT DISTINCT nx FROM Prop%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx']) 
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM Prop%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
        #initialise folders
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Square.txt','w')
        fid2=open(tmpdir + '/SquareProp.txt','w')
        
        #loop over columns
        s=0
        for nx in nxs:
            for ny in nys:
                req='SELECT p FROM Prop%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
                #print(req)
                self.fetch(req,['p'])
                ps=self.tmp.p
                square=False
                try:
                    if len(ps) == 4: 
                        #calculate_cm_orientation(self.tmp.p)
                        p=ps[0]            
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                        
                        Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                        Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                        
                        #Xvec=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                        #Yvec=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                        Xvec=[list(np.copy(self.Part.data.X)[:self.infos.tmax])]
                        Yvec=[list(np.copy(self.Part.data.Y)[:self.infos.tmax])]
                        
                        Phicm=np.zeros(self.infos.tmax)
                        for p in ps[1:]:
                            self.part(p)
                            self.Part.fields(["X","Y","t"])
                            Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                            Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                            Xvec+=[list(np.copy(self.Part.data.X)[:self.infos.tmax])]
                            Yvec+=[list(np.copy(self.Part.data.Y)[:self.infos.tmax])]
                        
                        Xcm=Xcm/4
                        Ycm=Ycm/4
                        Xvec=np.array(Xvec)
                        Yvec=np.array(Yvec)
                        
                        
                        Phicm=np.zeros(len(Xcm))
                        for ind in range(len(Phicm)):
                            
                            Phicm[ind] = np.arctan(      (np.average(Yvec[np.argwhere(Xvec[:,ind]>Xcm[ind]),ind]) -Ycm[ind])   /(np.average(Xvec[np.argwhere(Xvec[:,ind]>Xcm[ind]),ind]) -Xcm[ind]))
                        #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                        #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                        square=True
                                        
    
                    if len(ps) == 3: 
                        #calculate_cm_orientation(self.tmp.p)
                        p=ps[0]            
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                        
                        Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                        Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                        
                        Xvec=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                        Yvec=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                        
                        Phicm=np.zeros(self.infos.tmax)
                        for p in ps[1:]:
                            self.part(p)
                            self.Part.fields(["X","Y","t"])
                            Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                            Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                            Xvec+=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                            Yvec+=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                        
                        Xvec=np.array(Xvec)
                        Yvec=np.array(Yvec)
                        
                        dvec = [np.sqrt(    (Xvec[1]-Xvec[0])**2   +  (Yvec[1]-Yvec[0])**2),       np.sqrt(    (Xvec[2]-Xvec[0])**2   +  (Yvec[2]-Yvec[0])**2)  , np.sqrt(    (Xvec[2]-Xvec[1])**2   +  (Yvec[2]-Yvec[1])**2)   ]
                        Xmidvec = [(Xvec[1]+Xvec[0])/2,     (Xvec[2]+Xvec[0])/2,        (Xvec[2]+Xvec[1])/2     ]
                        Ymidvec = [(Yvec[1]+Yvec[0])/2,     (Yvec[2]+Yvec[0])/2,        (Yvec[2]+Yvec[1])/2     ]
                        otherindvec=[2,1,0]
                        
                        ind_diag = np.argmax(dvec)
                            
                        
                        Xcm = np.array([Xmidvec[np.argmax(dvec)]])
                        Ycm = np.array([Ymidvec[np.argmax(dvec)]])
                        
                        X4 = Xcm-(Xvec[otherindvec[ind_diag]]-Xcm)
                        Y4 = Ycm-(Yvec[otherindvec[ind_diag]]-Ycm)
                        
                        Xvec=np.concatenate((Xvec,X4),axis=0)
                        Yvec=np.concatenate((Yvec,Y4),axis=0)
                        
                        Phicm = np.arctan(      (np.average(Yvec[np.argwhere(Xvec>Xcm)]) -Ycm)   /(np.average(Xvec[np.argwhere(Xvec>Xcm)]) -Xcm))
                        
                        #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                        #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                        square=True
                    """
                    if len(ps) == 2: 
                        #calculate_cm_orientation(self.tmp.p)
                        p=ps[0]            
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                        
                        Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                        Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                        
                        Xvec=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                        Yvec=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                        
                        Phicm=np.zeros(self.infos.tmax)
                        for p in ps[1:]:
                            self.part(p)
                            self.Part.fields(["X","Y","t"])
                            Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                            Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                            Xvec+=list(np.copy(self.Part.data.X)[:self.infos.tmax])
                            Yvec+=list(np.copy(self.Part.data.Y)[:self.infos.tmax])
                        
                        Xvec=np.array(Xvec)
                        Yvec=np.array(Yvec)
                        
                        
                        dX= Xvec[1]-Xvec[0]
                        dY= Yvec[1]-Yvec[0]
                        
                        if nx ==nxs[0] and ny >nys[0]+0.1 and ny < nys[-1]+0.1:     #left edge
                            edgetype = 0
                            if dY<0:
                                dX=-dX
                                dY=-dY
                            Xvec2 = Xvec-dY
                            Yvec2 = Yvec+dX
                        elif nx ==nxs[-1] and ny >nys[0]+0.1 and ny < nys[-1]+0.1:        #right edge
                            edgetype = 1
                            if dY<0:
                                dX=-dX
                                dY=-dY
                            Xvec2 = Xvec+dY
                            Yvec2 = Yvec-dX                        
                        elif ny ==nys[0] and nx >nxs[0]+0.1 and nx < nxs[-1]+0.1:     #bottom edge
                            edgetype = 2
                            if dX<0:
                                dX=-dX
                                dY=-dY
                            Xvec2 = Xvec+dY
                            Yvec2 = Yvec-dX
                        elif ny ==nys[-1] and nx >nxs[0]+0.1 and nx < nxs[-1]+0.1:        #top edge
                            edgetype = 3
                            if dX<0:
                                dX=-dX
                                dY=-dY
                            Xvec2 = Xvec-dY
                            Yvec2 = Yvec+dX
                        else:
                            print('Error: can only track 2 holes but is not on edge')
                            ERRORlksemnrgv
                            
                        Xvec=np.concatenate((Xvec,Xvec2),axis=0)
                        Yvec=np.concatenate((Yvec,Yvec2),axis=0)
                            
                        Xcm = np.array([np.average(Xvec)])
                        Ycm = np.array([np.average(Yvec)])
                        
                        
                        Phicm = np.arctan(      (np.average(Yvec[np.argwhere(Xvec>Xcm)]) -Ycm)   /(np.average(Xvec[np.argwhere(Xvec>Xcm)]) -Xcm))
                        #Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                        #Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]
                        square=True
                    """
                            
                            
                        
                    """
                    if len(ps) == 2: 
                        
                        
                        
                        p=ps[0]            
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                        Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                        Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                        Phicm=np.zeros(self.infos.tmax)
                        p=ps[1]
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        
                        vecX=Xcm-np.copy(self.Part.data.X)[:self.infos.tmax]
                        vecY=Ycm-np.copy(self.Part.data.Y)[:self.infos.tmax]
                        
                        Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                        Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                        
                        offx=-np.abs(vecX)/2
                        offy=np.abs(vecY)/2
                        
                        Xcm=Xcm/2+offx
                        Ycm=Ycm/2+offy
                        
                        Phicm=np.arctan2(vecY,vecX)
                        Phicm+=-Phicm[0]#np.pi/4.-Phicm[0]            
                        square=True
                    """
                    #check
                    if ax1 is not None and square: 
                        ax1.scatter(Xcm,Ycm,color="r")
                     #Save for database
                    if square:
                        fid2.write('%d;%d;%d\n'% (s,nx,ny))
                        for ix in range(len(tcm)):
                            fid.write('%d;%d;%f;%f;%f\n'% (s,tcm[ix],Xcm[ix],Ycm[ix],Phicm[ix]))
                        s+=1
                
                except:
                    pass
        

        fid.close()
        fid2.close()
        
        req='DROP TABLE Square' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE SquareProp' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass        
        
        req='CREATE TABLE Square' + str(self.sub) + ' ( s int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL , Phi float NOT NULL , PRIMARY KEY ( t ,  s )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE SquareProp' + str(self.sub) + ' ( s int( 8  )  NOT  NULL ,  nx int( 8  )  NOT  NULL , ny int( 8  )  NOT  NULL , PRIMARY KEY ( s , nx, ny )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Square.txt' INTO TABLE Square" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/SquareProp.txt' INTO TABLE SquareProp" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)        
        
        
        
        
    def writetxtSQUARE_backup(self,ax1=None):          
        #retrieve data
        req='SELECT DISTINCT nx FROM Prop%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx']) 
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM Prop%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
        #initialise folders
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Square.txt','w')
        fid2=open(tmpdir + '/SquareProp.txt','w')
        
        #loop over columns
        s=0
        for nx in nxs:
            for ny in nys:
                #retrieve particle label corresponding to target nx and ny
                req='SELECT p FROM Prop%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
                self.fetch(req,['p'])
                ps=self.tmp.p
                
                #calculate square based on positions of its points
                square=False
                if len(ps) == 4:# if square is made of 4 points
                    #first point
                    p=ps[0]            
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    Phicm=np.zeros(self.infos.tmax)
                    #other points
                    for p in ps[1:]:
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                        Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                    #Calculate center of mass and orientation of square
                    Xcm=Xcm/4
                    Ycm=Ycm/4
                    Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                    Phicm+=np.pi/4.-Phicm[0]
                    square=True
                
                if len(ps) == 2:# if square is made of 2 points  
                    #first point
                    p=ps[0]            
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    Phicm=np.zeros(self.infos.tmax)
                    #second point
                    p=ps[1]
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    #calculate center of mass and orientation
                    vecX=Xcm-np.copy(self.Part.data.X)[:self.infos.tmax]
                    vecY=Ycm-np.copy(self.Part.data.Y)[:self.infos.tmax]
                    
                    Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                    Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                    
                    #!!!!! this depends on where the 2 points are located
                    offx=-np.abs(vecX)/2
                    offy=np.abs(vecY)/2
                    
                    Xcm=Xcm/2+offx
                    Ycm=Ycm/2+offy
                    
                    Phicm=np.arctan2(vecY,vecX)
                    Phicm+=np.pi/4.-Phicm[0]            
                    square=True
                
                #visual check
                if ax1 is not None and square: 
                    ax1.scatter(Xcm,Ycm,color="r")
                #Save for database
                if square:
                    fid2.write('%d;%d;%d\n'% (s,nx,ny))
                    for ix in range(len(tcm)):
                        fid.write('%d;%d;%f;%f;%f\n'% (s,tcm[ix],Xcm[ix],Ycm[ix],Phicm[ix]))
                    s+=1
        fid.close()
        fid2.close()
        
        #reinitialize tables in DB
        #remove
        req='DROP TABLE Square' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE SquareProp' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass        
        #initialize table format
        req='CREATE TABLE Square' + str(self.sub) + ' ( s int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL , Phi float NOT NULL , PRIMARY KEY ( t ,  s )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE SquareProp' + str(self.sub) + ' ( s int( 8  )  NOT  NULL ,  nx int( 8  )  NOT  NULL , ny int( 8  )  NOT  NULL , PRIMARY KEY ( s , nx, ny )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        
        #load data into DB
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Square.txt' INTO TABLE Square" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/SquareProp.txt' INTO TABLE SquareProp" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)        
        
        
        
    def writetxtELLIPSE(self,ax1=None):
        #select nx and ny of ellipses
        req='SELECT DISTINCT nx FROM SquareProp%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx'])
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM SquareProp%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
        
        #loop overs ellipses
        ee=[]
        for nx in nxs:
            for ny in nys:
                #Select diagonally arranged squares
                
                
                #m = middle
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx,ny)
                self.fetch(req,['p','nx','ny'])
                pbm=np.copy(self.tmp.p);nxbm=np.copy(self.tmp.nx);nybm=np.copy(self.tmp.ny)
                
                #l = left
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx-1,ny)
                self.fetch(req,['p','nx','ny'])
                pbl=np.copy(self.tmp.p);nxbl=np.copy(self.tmp.nx);nybl=np.copy(self.tmp.ny)
                
                #r = right
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx+1,ny)
                self.fetch(req,['p','nx','ny'])
                pbr=np.copy(self.tmp.p);nxbr=np.copy(self.tmp.nx);nybr=np.copy(self.tmp.ny)
                
                #b = bottom
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx,ny-1)
                self.fetch(req,['p','nx','ny'])
                pbb=np.copy(self.tmp.p);nxbb=np.copy(self.tmp.nx);nybb=np.copy(self.tmp.ny)
                
                #t = top
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx,ny+1)
                self.fetch(req,['p','nx','ny'])
                pbt=np.copy(self.tmp.p);nxbt=np.copy(self.tmp.nx);nybt=np.copy(self.tmp.ny)
                

                print(pbm,pbl,pbr,pbb,pbt)       
                #either left vs right or top  vs bottom
                if len(pbl)*len(pbr)==1 and len(pbm)==0:
                    ee+=[["leftright",pbl[0],pbr[0],nx,ny]]
                elif len(pbb)*len(pbt)==1 and len(pbm)==0:
                    ee+=[["bottomtop",pbb[0],pbt[0],nx,ny]]
                    

        #Define ellipses from the squares
        patches = []
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Ellipse.txt','w')
        fid2=open(tmpdir + '/EllipseProp.txt','w')
        
        for ie,e in zip(range(len(ee)),ee):
        #req='SELECT s,X,Y FROM Square%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nxp,ny)        
        #self.fetch(req,['p','nx','ny'])
        #pbr=np.copy(self.tmp.p);nxbr=np.copy(self.tmp.nx);nybr=np.copy(self.tmp.ny)
            req='SELECT s,X,Y,Phi,t FROM Square%s WHERE (s=%s) ORDER BY t;' % (self.sub,e[1])
            self.fetch(req,['s','X','Y','Phi','t'])
            
            t1=np.copy(self.tmp.t);
            X1=np.copy(self.tmp.X);Y1=np.copy(self.tmp.Y);Phi1=np.copy(self.tmp.Phi)
            
            req='SELECT s,X,Y,Phi,t FROM Square%s WHERE (s=%s) ORDER BY t;' % (self.sub,e[2])
            self.fetch(req,['s','X','Y','Phi','t'])
            
            t2=np.copy(self.tmp.t);
            X2=np.copy(self.tmp.X);Y2=np.copy(self.tmp.Y);Phi2=np.copy(self.tmp.Phi)
            
            
            lt = np.min([len(t1),len(t2)])
            
            
            t1=t1[:lt];     X1=X1[:lt];     Y1=Y1[:lt];     Phi1=Phi1[:lt];
            t2=t2[:lt];     X2=X2[:lt];     Y2=Y2[:lt];     Phi2=Phi2[:lt];
            
            
            
            ellipse_Xcm=(X1+X2)/2
            ellipse_Ycm=(Y1+Y2)/2
            boxX1=[X1+np.sqrt(2)/2*np.cos(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]     #boxX1=[X1+0.5*np.cos(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
            boxY1=[Y1+np.sqrt(2)/2*np.sin(Phi1-np.pi/4+i*np.pi/2) for i in range(5)] # boxY1=[Y1+0.5*np.sin(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
            
            boxX2=[X2+np.sqrt(2)/2*np.cos(Phi2-np.pi/4+i*np.pi/2) for i in range(5)]#!!!         #boxX2=[X2+0.5*np.cos(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]#!!!
            boxY2=[Y2+np.sqrt(2)/2*np.sin(Phi2-np.pi/4+i*np.pi/2) for i in range(5)]         #boxY2=[Y2+0.5*np.sin(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
        
            
            if e[0] == "leftright":
                ellipse_vecXH=boxX2[3]-boxX1[1]
                ellipse_vecYH=boxY2[3]-boxY1[1]
                
                ellipse_vecXV=boxX2[2]-boxX1[0]
                ellipse_vecYV=boxY2[2]-boxY1[0]
                
                
            elif e[0] == "bottomtop":
                ellipse_vecXH=boxX2[0]-boxX1[2]
                ellipse_vecYH=boxY2[0]-boxY1[2]
                
                ellipse_vecXV=-boxX2[3]+boxX1[1]
                ellipse_vecYV=-boxY2[3]+boxY1[1]
            
            a=np.sqrt(ellipse_vecXH**2+ellipse_vecYH**2)
            b=np.sqrt(ellipse_vecXV**2+ellipse_vecYV**2)
            angle=np.arctan2(ellipse_vecYH,ellipse_vecXH)
            
            
            
            fid2.write('%d;%d;%d;%d;%d;%s\n'% (ie,e[1],e[2],e[3],e[4],e[0]))
            for ix in range(len(t1)):
                fid.write('%d;%d;%f;%f;%f;%f;%f\n'% (ie,t1[ix],ellipse_Xcm[ix],ellipse_Ycm[ix],a[ix],b[ix],angle[ix]))
            
            if ax1 is not None:
                t=0
                ax1.plot( [boxX1[i][t] for i in range(5)],[boxY1[i][t] for i in range(5)],c="k",lw=1)
                ax1.plot( [boxX2[i][t] for i in range(5)],[boxY2[i][t] for i in range(5)],c="k",lw=1)    
                
                
                el=Ellipse(xy=[ellipse_Xcm[t],ellipse_Ycm[t]],
                                   width=a[t],height=b[t],angle=angle[t])
                
                patches.append(el)
        if ax1 is not None:
            p = PatchCollection(patches, alpha=1,facecolor="None",edgecolor="g",linewidth=2)
            ax1.add_collection(p)
        fid.close()
        fid2.close()
        
        req='DROP TABLE Ellipse' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE EllipseProp' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass        
        
        req='CREATE TABLE Ellipse' + str(self.sub) + ' ( e int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL, a float NOT  NULL , b float NOT NULL , angle float NOT NULL , PRIMARY KEY ( t ,  e )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE EllipseProp' + str(self.sub) + ' ( e int( 8  )  NOT  NULL ,  p1 int( 8  )  NOT  NULL , p2 int( 8  )  NOT  NULL,  nx int( 8  )  NOT  NULL , ny int( 8  )  NOT  NULL ,typelink char( 5 ) NOT  NULL, PRIMARY KEY ( e , nx, ny )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Ellipse.txt' INTO TABLE Ellipse" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/EllipseProp.txt' INTO TABLE EllipseProp" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        
             
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    def writetxtELLIPSE_backup(self,ax1=None):
        #select nx and ny of ellipses
        req='SELECT DISTINCT nx FROM SquareProp%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx'])
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM SquareProp%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
        
        #loop overs ellipses
        ee=[]
        for nx,nxp in zip(nxs[:-1],nxs[1:]):
            for ny,nyp in zip(nys[:-1],nys[1:]):
                #Select diagonally arranged squares
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx,ny)
                self.fetch(req,['p','nx','ny'])
                pbl=np.copy(self.tmp.p);nxbl=np.copy(self.tmp.nx);nybl=np.copy(self.tmp.ny)
                
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nxp,ny)        
                self.fetch(req,['p','nx','ny'])
                pbr=np.copy(self.tmp.p);nxbr=np.copy(self.tmp.nx);nybr=np.copy(self.tmp.ny)
                
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nx,nyp)        
                self.fetch(req,['p','nx','ny'])
                ptl=np.copy(self.tmp.p);nxtl=np.copy(self.tmp.nx);nytl=np.copy(self.tmp.ny)
                
                req='SELECT s,nx,ny FROM SquareProp%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nxp,nyp)        
                self.fetch(req,['p','nx','ny'])
                ptr=np.copy(self.tmp.p);nxtr=np.copy(self.tmp.nx);nytr=np.copy(self.tmp.ny)
        
                print(pbl,pbr,ptl,ptr)       
                #either bottomleft vs topright or top left vs bottomright
                if len(pbl)*len(ptr)==1:
                    ee+=[["diag1",pbl[0],ptr[0],nx,ny]]
                if len(ptl)*len(pbr)==1:
                    ee+=[["diag2",ptl[0],pbr[0],nx,ny]]
        
        #Define ellipses from the squares
        patches = []
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Ellipse.txt','w')
        fid2=open(tmpdir + '/EllipseProp.txt','w')
        
        for ie,e in zip(range(len(ee)),ee):
        #req='SELECT s,X,Y FROM Square%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nxp,ny)        
        #self.fetch(req,['p','nx','ny'])
        #pbr=np.copy(self.tmp.p);nxbr=np.copy(self.tmp.nx);nybr=np.copy(self.tmp.ny)
            req='SELECT s,X,Y,Phi,t FROM Square%s WHERE (s=%s) ORDER BY t;' % (self.sub,e[1])
            self.fetch(req,['s','X','Y','Phi','t'])
            
            t1=np.copy(self.tmp.t);
            X1=np.copy(self.tmp.X);Y1=np.copy(self.tmp.Y);Phi1=np.copy(self.tmp.Phi)
            
            req='SELECT s,X,Y,Phi,t FROM Square%s WHERE (s=%s) ORDER BY t;' % (self.sub,e[2])
            self.fetch(req,['s','X','Y','Phi','t'])
            
            t2=np.copy(self.tmp.t);
            X2=np.copy(self.tmp.X);Y2=np.copy(self.tmp.Y);Phi2=np.copy(self.tmp.Phi)
            
            ellipse_Xcm=(X1+X2)/2
            ellipse_Ycm=(Y1+Y2)/2
            boxX1=[X1+0.5*np.cos(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
            boxY1=[Y1+0.5*np.sin(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
            
            #!!!CHECK !!!
            boxX2=[X2+0.5*np.cos(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]#!!!
            #!!!CHECK !!!
            boxY2=[Y2+0.5*np.sin(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
        
            if e[0] == "diag1":
                ellipse_vecXH=boxX2[3]-boxX1[1]
                ellipse_vecYH=boxY2[3]-boxY1[1]
                
                ellipse_vecXV=boxX2[2]-boxX1[0]
                ellipse_vecYV=boxY2[2]-boxY1[0]
            elif e[0] == "diag2":
                ellipse_vecXH=boxX2[1]-boxX1[3]
                ellipse_vecYH=boxY2[1]-boxY1[3]
                
                ellipse_vecXV=-boxX2[2]+boxX1[0]
                ellipse_vecYV=-boxY2[2]+boxY1[0]
            
            a=np.sqrt(ellipse_vecXH**2+ellipse_vecYH**2)
            b=np.sqrt(ellipse_vecXV**2+ellipse_vecYV**2)
            angle=np.arctan2(ellipse_vecYH,ellipse_vecXH)
            
            fid2.write('%d;%d;%d;%d;%d;%s\n'% (ie,e[1],e[2],e[3],e[4],e[0]))
            for ix in range(len(t1)):
                fid.write('%d;%d;%f;%f;%f;%f;%f\n'% (ie,t1[ix],ellipse_Xcm[ix],ellipse_Ycm[ix],a[ix],b[ix],angle[ix]))
            
            if ax1 is not None:
                t=0
                ax1.plot( [boxX1[i][t] for i in range(5)],[boxY1[i][t] for i in range(5)],c="k",lw=1)
                ax1.plot( [boxX2[i][t] for i in range(5)],[boxY2[i][t] for i in range(5)],c="k",lw=1)    
                el=Ellipse(xy=[ellipse_Xcm[t],ellipse_Ycm[t]],
                                   width=a[t],height=b[t],angle=angle[t])
                patches.append(el)
        if ax1 is not None:
            p = PatchCollection(patches, alpha=1,facecolor="None",edgecolor="g",linewidth=2)
            ax1.add_collection(p)
        fid.close()
        fid2.close()
        
        req='DROP TABLE Ellipse' + str(self.sub)
        try:
            self.DB.cursor.execute(req)
        except: pass
        req='DROP TABLE EllipseProp' + str(self.sub)
        try:        
            self.DB.cursor.execute(req)
        except: pass        
        
        req='CREATE TABLE Ellipse' + str(self.sub) + ' ( e int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL, a float NOT  NULL , b float NOT NULL , angle float NOT NULL , PRIMARY KEY ( t ,  e )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        req='CREATE TABLE EllipseProp' + str(self.sub) + ' ( e int( 8  )  NOT  NULL ,  p1 int( 8  )  NOT  NULL , p2 int( 8  )  NOT  NULL,  nx int( 8  )  NOT  NULL , ny int( 8  )  NOT  NULL ,typelink char( 5 ) NOT  NULL, PRIMARY KEY ( e , nx, ny )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
        self.DB.cursor.execute(req)
        
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/Ellipse.txt' INTO TABLE Ellipse" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        req="LOAD DATA LOCAL INFILE '" + tmpdir + "/EllipseProp.txt' INTO TABLE EllipseProp" +  str(self.sub) + " FIELDS TERMINATED BY ';' LINES TERMINATED BY '\n'"
        self.DB.cursor.execute(req)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    def createTables(self):
        req='CREATE TABLE Data' + str(self.sub) + ' ( p int( 8  )  NOT  NULL , t mediumint( 6  )  NOT  NULL , X float NOT  NULL , Y float NOT NULL , area float NOT NULL, PRIMARY KEY (  t ,  p  )    ) ENGINE  =  MyISAM  DEFAULT CHARSET  = utf8;'
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
        print(req)
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
        self.path=path_sys[_platform]+'Experiments/'
        self.name=name
        self.sub=sub
        self.DB=Database.DB()
        self.DB.connect(self.name)
        self.tmp=Tmp()
        
        self.infos=Infos()
        self.infos.file=self.path + self.name + '/Files/infos' + str(self.sub) + '.p'
        
        self.imagedata=Imagedata()        
        self.imagedata.file=self.path + self.name + '/Files/imagedata' + str(self.sub) + '.p'

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
        #self.plotsquares(ax1,t,p=[pi1],color="b")
        ax1.scatter(X1[pi1],Y1[pi1],color="r")

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
        #self.plotsquares(ax1,tb,color="b")
        self.Time.fields(['t','Y','X','p'])
        ax1.scatter(self.Time.data.X,self.Time.data.Y,color="b")
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
        #self.plotsquares(ax1,t,p=[pi1],color="r",alpha=1)
        ax1.scatter(self.Time.data.X,self.Time.data.Y,color="r")
        pi2,ix2=Tools.pick(p2,X2,Y2)
        ax1.scatter()
        #self.plotsquares(ax1,tb,p=[pi2],color="r",alpha=1)
        
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
        
    def assign_nxny_auto(self,ax1,offx=0,offy=0):
        X1=self.Time.data.X;Y1=self.Time.data.Y;p1=self.Time.data.p;
        
        for i in range(self.infos.N):
            x1=i+offx
            x2=x1+1
            y1=i+offy
            y2=y1+1
            ax1.axvline(x=x1,lw=1)
            ax1.axhline(y=y1,lw=1)
            
            condX1=np.double((X1>=x1)==True)
            condX2=np.double((X1<=x2)==True)
                                    
            ixX=find(condX1+condX2==2)
            ixX=np.array(ixX,int)
            poutX=p1[ixX]        
            
            if len(poutX)>0.5:
                self.updatenx(poutX,i+1)
            
            condY1=np.double((Y1>=y1)==True)
            condY2=np.double((Y1<=y2)==True)
            
            ixY=find(condY1+condY2==2)
            ixY=np.array(ixY,int)
            poutY=p1[ixY]        
            if len(poutY)>0.5:
                self.updateny(poutY,i+1)
        ax1.axvline(x=x2,lw=1)
        ax1.axhline(y=y2,lw=1)
        ax1.scatter(self.Time.data.X,self.Time.data.Y)
        x=plt.ginput()
