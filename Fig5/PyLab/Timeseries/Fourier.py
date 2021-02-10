# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 15:23:52 2014
Home made fourier serie decomposition

@author: coulais
"""
import numpy as np

class Fseries:
    
    def __init__(self,y,N):
        self.y=y
        self.N=N
        self.Np=len(self.y)
        self.x=np.linspace(0.,1.,self.Np)
        self.can=np.zeros(N)
        self.cbn=np.zeros(N)
        self.decomposition()
    
    def an(self,H):
        if H==0:
            c = self.y*np.cos(0.0*self.x)/2.#*Np
        else:
            c = (-1)**0*self.y*np.cos(2*H*np.pi*self.x)#*Np
        self.can[H]=np.trapz(c,dx=1./float(self.Np))*2.
        
    def bn(self,H):
        c = (-1)**0*self.y*np.sin(2*H*np.pi*self.x)#*Np
        self.cbn[H]=np.trapz(c,dx=1./float(self.Np))*2.    
    
    def decomposition(self):
        self.s=np.zeros(len(self.y))
        self.sins=[np.zeros(len(self.y))]*self.Np
        self.coss=[np.zeros(len(self.y))]*self.Np

        for H in range(self.N):
            self.an(H)
            self.bn(H)
            self.s+=self.can[H]*np.cos(2*H*np.pi*self.x)+self.cbn[H]*np.sin(2*H*np.pi*self.x)
            self.coss[H]=self.can[H]*np.cos(2*H*np.pi*self.x)
            self.sins[H]=self.cbn[H]*np.sin(2*H*np.pi*self.x)
