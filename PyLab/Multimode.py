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
import imutils
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
path_sys={"darwin":"/Users/coulais/science/Shared/Metacombinatorial/",
          "win32":"TBD",
          "linux":"TBD"}


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
        #to comment if no mysql server        
        self.DB=Database.DB()
        self.DB.connect(self.name)
        self.tmp=Tmp()
        #to comment if no mysql server
        
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
        for ff in os.listdir(path+'Frames'):
            if ".tiff" in ff:
                num+=[int(ff[6:-5])]
                #print(num)
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
        for jj in pics:
            file_in= "frame_%04d.tiff" % (jj)
            file_out= "frame_%04d.jpg" % (jj)
            #Load image
            print(path +'Frames/'+ file_in)
            originalImage = cv2.imread(path + file_in)
            #img = cv2.imread(path+'Frames/' + file_in,0)
            try:
            #if originalImage is not None: #Does the image exist?
                #///////////////////// INITIAL STEP ///////////////////////
                ######## Load ####
                img = cv2.imread(path+'Frames/' + file_in,0)#load image as grayscale image
                back = cv2.imread(path+'Background/' + file_in,0)
                CLONE = np.copy(img)      # to draw check picture  
                CLONE3 = np.copy(originalImage)
                #CLONE2 = cv2.cvtColor(CLONE,cv2.COLOR_RGB2GRAY) 
                ######## Rotate and Crop ####
                x0,deltax,y0,deltay=self.infos.BBox
                imgcrop = np.copy(img[y0:y0+deltay, x0:x0+deltax])#Crop image)
                backcrop = np.copy(back[y0:y0+deltay, x0:x0+deltax])#Crop image)
                del originalImage
                #imgray = imgcrop.astype('float32') 
                #backgray = backcrop.astype('float32') 
                imgray = imgcrop.astype('int16') 
                backgray = backcrop.astype('int16') 
                img2=(imgray-backgray)/2+125
                img2 = img2.astype('uint8')                 
                if display is True:
                    ax1.imshow(img2)
                (tmp,ThresholdedImage)= cv2.threshold(img2,self.infos.ths,255,cv2.THRESH_BINARY_INV)
                #(tmp,ThresholdedImage)= cv2.threshold(imgcrop,255-self.infos.ths2,255-self.infos.ths1,cv2.THRESH_BINARY_INV)
                del tmp
                
                #ax2.imshow(ThresholdedImage)
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
                        if display is True:                        
                            ax2.plot(x,y,marker=".",c="w",markersize=12,markeredgewidth=6)
                        cv2.drawContours(imgcrop,[cc],-1,(125,255,125),-1)
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
        file_in= "frame_%04d.tiff" % (jj)
        img = cv2.imread(path+'Frames/' + file_in,0)
        x0,deltax,y0,deltay=self.infos.BBox
        imgcrop = np.copy(img[y0:y0+deltay, x0:x0+deltax])#Crop 
        ax1.imshow(imgcrop,cmap="gray")                        
        scaling=1/self.infos.pix2mm*self.infos.square
        sizes=0.8
        offset=[0,0]#[self.infos.BBox[0],self.infos.BBox[2]]
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
            
            #el.set_alpha(0.5)
            ialpha=np.abs(Omega[e])*2
            if ialpha >=1: ialpha=1.0
            #ialpha=1
            el.set_alpha(ialpha)
            icol=(Omega[e]*(-1)**(nxref[ix][0]+nyref[ix][0])+1)/2.
            el.set_facecolor(plt.cm.bwr(icol))#plt.cm.PuOr(icol))
            el.set_edgecolor("None")
            ax1.add_artist(el)
            
        '''
        req='SELECT DISTINCT s FROM SquareProp%s ;' % (self.sub)
        self.fetch(req,['s'])
        ss=np.copy(self.tmp.s)
        
        self.Time.fields(['t','Y','X','p'])
        ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[0])
        
        for ie,e in zip(range(len(ss)),ss):
        #req='SELECT s,X,Y FROM Square%s WHERE ((nx=%s AND ny=%s));' % (self.sub,nxp,ny)        
        #self.fetch(req,['p','nx','ny'])
        #pbr=np.copy(self.tmp.p);nxbr=np.copy(self.tmp.nx);nybr=np.copy(self.tmp.ny)
            req='SELECT s,X,Y,Phi,t FROM Square%s WHERE (t=%s AND s=%s) ;' % (self.sub,self.Time.t,e)
            #print(req)
            self.fetch(req,['s','X','Y','Phi','t'])
            
            t1=np.copy(self.tmp.t);
            X1=np.copy(self.tmp.X);Y1=np.copy(self.tmp.Y);Phi1=np.copy(self.tmp.Phi)
        
            ax1.scatter(X1*scaling+offset[0],Y1*scaling+offset[0],c="r")
    '''
        
            
    def extract_ellipses_orderparam(self,t,nx,ny):
        #retrieve ellipses properties
        self.time(t)
        #self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT e,nx,ny FROM EllipseProp%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
        self.fetch(req,['e','nx','ny'])
        #print(req)
        eref=np.copy(self.tmp.e);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        #for the time step t
        try:
            req='SELECT e,X,Y,a,b,angle,t FROM Ellipse%s WHERE (t=%s AND e=%s);' % (self.sub,self.Time.t,eref[0])
            #print(req)
            self.fetch(req,['e','X','Y','a','b','angle','t'])            
            t=np.copy(self.tmp.t);elab=np.copy(self.tmp.e);
            X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);a=np.copy(self.tmp.a);b=np.copy(self.tmp.b);angle=np.copy(self.tmp.angle)
            
            #plot ellipses
            self.Time.data.ny=nyref
            self.Time.data.nx=nxref
            self.Time.data.e=elab
            self.Time.data.X=X
            self.Time.data.Y=Y
            self.Time.data.a=a
            self.Time.data.b=b
            self.Time.data.angle=angle
            self.Time.data.Omega=(1-a/b)#*np.cos(2*angle)
        except:
            self.Time.data.ny=np.NaN
            self.Time.data.nx=np.NaN
            self.Time.data.e=np.NaN
            self.Time.data.X=np.NaN
            self.Time.data.Y=np.NaN
            self.Time.data.a=np.NaN
            self.Time.data.b=np.NaN
            self.Time.data.angle=np.NaN
            self.Time.data.Omega=np.NaN
    def extract_strain(self,t,nx,ny):
        #retrieve ellipses properties
        self.time(t)
        #self.Time.fields(['t','Y','X','p'])
        #ax1.scatter(self.Time.data.X*scaling+offset[0],self.Time.data.Y*scaling+offset[1])
        req='SELECT s,nx,ny FROM SquareProp%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
        self.fetch(req,['s','nx','ny'])
        #print(req)
        sref=np.copy(self.tmp.s);nxref=np.copy(self.tmp.nx);nyref=np.copy(self.tmp.ny);
        
        #for the time step t
        try:
            req='SELECT s,X,Y,Phi,t FROM Square%s WHERE (t=%s AND s=%s);' % (self.sub,self.Time.t,sref[0])
            #print(req)
            self.fetch(req,['s','X','Y','Phi','t'])            
            t=np.copy(self.tmp.t);slab=np.copy(self.tmp.s);
            X=np.copy(self.tmp.X);Y=np.copy(self.tmp.Y);Phi=np.copy(self.tmp.Phi)
            
            #plot ellipses
            self.Time.data.ny=nyref
            self.Time.data.nx=nxref
            self.Time.data.s=slab
            self.Time.data.X=X
            self.Time.data.Y=Y
            self.Time.data.Phi=Phi

        except:
            self.Time.data.ny=np.NaN
            self.Time.data.nx=np.NaN
            self.Time.data.s=np.NaN
            self.Time.data.X=np.NaN
            self.Time.data.Y=np.NaN
            self.Time.data.Phi=np.NaN


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
        req='SELECT DISTINCT nx FROM Prop%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx'])
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM Prop%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
        tmpdir=self.path + self.name + '/tmpDB/' + str(self.sub) 
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)        
        fid=open(tmpdir + '/Square.txt','w')
        fid2=open(tmpdir + '/SquareProp.txt','w')
        
        s=0
        for nx in nxs:
            for ny in nys:
                req='SELECT p FROM Prop%s WHERE (nx=%s AND ny=%s);' % (self.sub,nx,ny)
                #print(req)
                self.fetch(req,['p'])
                ps=self.tmp.p
                square=False
                if len(ps) == 4: 
                    #calculate_cm_orientation(self.tmp.p)
                    p=ps[0]            
                    self.part(p)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    Xcm=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Ycm=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    Phicm=np.zeros(self.infos.tmax)
                    for p in ps[1:]:
                        self.part(p)
                        self.Part.fields(["X","Y","t"])
                        Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                        Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                    Xcm=Xcm/4
                    Ycm=Ycm/4
                    Phicm=np.arctan2(self.Part.data.Y[:self.infos.tmax]-Ycm,self.Part.data.X[:self.infos.tmax]-Xcm)
                    Phicm+=np.pi/4.-Phicm[0]
                    square=True
                if len(ps) == 2: 
                    p1=ps[0]            
                    self.part(p1)
                    self.Part.fields(["X","Y","t"])
                    tcm=np.copy(self.Part.data.t)[:self.infos.tmax];
                    X1=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Y1=np.copy(self.Part.data.Y)[:self.infos.tmax];
                    Phicm=np.zeros(self.infos.tmax)
                    p2=ps[1]
                    self.part(p2)
                    self.Part.fields(["X","Y","t"])
                    X2=np.copy(self.Part.data.X)[:self.infos.tmax];
                    Y2=np.copy(self.Part.data.Y)[:self.infos.tmax];

                    vecX=X1-X2
                    vecY=Y1-Y2
                    #print(vecX,vecY)
                    
                    #Xcm+=np.copy(self.Part.data.X)[:self.infos.tmax]
                    #Ycm+=np.copy(self.Part.data.Y)[:self.infos.tmax]
                    
                    offx=-np.abs(vecX)/2
                    offy=-np.abs(vecY)/2
                    
                    Xcm=(X1+X2)/2+offx
                    Ycm=(Y1+Y2)/2+offy
                    
                    Phicm=np.arctan2(vecY,vecX)
                    Phicm+=np.pi/4.-Phicm[0]            
                    square=True
                #check
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
        
    def writetxtELLIPSE(self,ax1=None):   
        req='SELECT DISTINCT nx FROM SquareProp%s WHERE nx!=0 ;' % (self.sub)
        self.fetch(req,['nx'])
        nxs=np.sort(np.copy(self.tmp.nx))
        req='SELECT DISTINCT ny FROM SquareProp%s WHERE ny!=0 ;' % (self.sub)
        self.fetch(req,['ny'])
        nys=np.sort(np.copy(self.tmp.ny))
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
            
            boxX2=[X2+0.5*np.cos(Phi1-np.pi/4+i*np.pi/2) for i in range(5)]
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
        self.plotsquares(ax1,tb,p=[pi2],color="r",alpha=1)
        
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
            self.updatenx(poutX,i+1)
            
            condY1=np.double((Y1>=y1)==True)
            condY2=np.double((Y1<=y2)==True)
            
            ixY=find(condY1+condY2==2)
            ixY=np.array(ixY,int)
            poutY=p1[ixY]        
            self.updateny(poutY,i+1)
        ax1.axvline(x=x2,lw=1)
        ax1.axhline(y=y2,lw=1)
        ax1.scatter(self.Time.data.X,self.Time.data.Y)
        x=plt.ginput()
