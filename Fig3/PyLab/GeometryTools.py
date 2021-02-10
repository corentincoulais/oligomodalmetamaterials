# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 11:04:27 2016

@author: coulais
"""


import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import scipy.spatial

mpl.rcParams['axes.labelsize'] = 28
mpl.rcParams['xtick.labelsize'] =  28
mpl.rcParams['ytick.labelsize'] = 28
mpl.rcParams['lines.linewidth']=3
mpl.rcParams['lines.markersize']=10
mpl.rcParams['axes.linewidth']= 2.0
#mpl.rcParams['axes.formatter.limits']=[ -2, 2] )
#mpl.rcParams['legend.fontsize']=14
mpl.rcParams['font.family']='serif'
mpl.r/Volumes/data/AMOLF/groups/hecke-group/Corentin/rawdata/PyLocal/GeometryTools.pycParams['legend.fontsize']=22;mpl.rcParams['lines.linewidth']=2

def find_neighbors(pindex, triang):
    return triang.vertex_neighbor_vertices[1][triang.vertex_neighbor_vertices[0][pindex]:triang.vertex_neighbor_vertices[0][pindex+1]]

def angular_correlation(xy,ax=None,visualcheck=False):
    ###############################################################################
    ##### 1rst approach: correlation between neighbors distance and angle #########
    ###############################################################################
    ############## Extract Neighbors via voronoi tesselation #############
    vor = Voronoi(xy)
    if visualcheck == True:    
        voronoi_plot_2d(vor)
    ######## Calculate simple correlation over neighbors ########
    DR=np.array([])
    DT=np.array([])
    for p1,p2 in vor.ridge_points:#loop over neighbors
        dxy=xy[p1]-xy[p2]
        theta=np.arctan2(dxy[1],dxy[0])
        dr=np.sqrt(dxy[0]**2+dxy[1]**2)
        DR=np.hstack((DR,dr))
        DT=np.hstack((DT,theta))
    if ax is not None:
        #### Show scatter plot of neighbor  angle vs. distance ###
        ax.plot(DT, DR,linestyle="none",marker="d")

def angular_order_parameter(x_list,y_list,visualcheck=False,factor=8,limn=4):
    ###############################################################################
    ######## 2nd approach: calculate orientational order parameter#################
    ###############################################################################
    
    ######## Extract Neighbors via Delaunay triangulation ########
    tri = scipy.spatial.Delaunay(np.array([[x,y] for x,y in zip(x_list, y_list)]))
    
    ######## Calculate orientational order parameter ########
    DR=np.array([])
    DT=np.array([])
    PSIlist=[[]]*len(x_list)
    PSI=0
    PSII=0
    N=0
    
    for pindex in range(len(x_list)):
        neighbor_indices = find_neighbors(pindex,tri)    
        if visualcheck == True:
            plt.plot(x_list[pindex], y_list[pindex], marker='d',color="DeepPink")
            plt.plot([x_list[i] for i in neighbor_indices],
                       [y_list[i] for i in neighbor_indices], 'ro')    
            plt.show()
            
        dx=x_list[pindex]-x_list[neighbor_indices]
        dy=y_list[pindex]-y_list[neighbor_indices]
        theta=np.arctan2(dy,dx)
        dr=np.sqrt(dx**2+dy**2)
        DR=np.hstack((DR,dr))
        DT=np.hstack((DT,theta))
        #print len(theta)
        PSIlist[pindex]=np.sqrt(np.real(np.sum(np.exp(factor*1j*theta))/len(theta))**2+np.imag(np.sum(np.exp(factor*1j*theta))/len(theta))**2)
        if len(neighbor_indices) >= limn:
            PSI+=np.real(np.sum(np.exp(factor*1j*theta))/len(theta))
            PSII+=np.imag(np.sum(np.exp(factor*1j*theta))/len(theta))
            N+=1
    #N=len(x_list)
    PSI=PSI/N
    PSII=PSII/N
    return PSI,PSII,N,PSIlist

def angular_order_parameter2(x_list,y_list,visualcheck=False,factor=8,limn=4,ax=None):
    ###############################################################################
    ######## 2nd approach: calculate orientational order parameter#################
    ###############################################################################
    
    ######## Extract Neighbors via Delaunay triangulation ########
    #tri = scipy.spatial.Delaunay(np.array([[x,y] for x,y in zip(x_list, y_list)]))
    xy=np.vstack((np.array(x_list),np.array(y_list))).transpose()
    vor = Voronoi(xy,furthest_site=False)#,qhull_options="Qc")   
    neighbors=[[]]*len(x_list)
    for pindex in range(len(x_list)):
        neighbors[pindex]=vor.ridge_points[np.where(vor.ridge_points[:,0]==pindex),1].tolist()[0]+vor.ridge_points[np.where(vor.ridge_points[:,1]==pindex),0].tolist()[0]

    ######## Calculate orientational order parameter ########
    DR=np.array([])
    DT=np.array([])
    PSIlist=[[]]*len(x_list)
    PSI=0
    PSII=0
    N=0
    
    for pindex in range(len(neighbors)):
        #neighbor_indices = find_neighbors(pindex,tri)  
        neighbor_indices = neighbors[pindex]
        if False:#visualcheck == True:
            plt.plot(x_list[pindex], y_list[pindex], marker='d',color="DeepPink")
            plt.plot([x_list[i] for i in neighbor_indices],
                       [y_list[i] for i in neighbor_indices], 'ro')    
            plt.show()
            
        dx=x_list[pindex]-x_list[neighbor_indices]
        dy=y_list[pindex]-y_list[neighbor_indices]
        theta=np.arctan2(dy,dx)
        dr=np.sqrt(dx**2+dy**2)
        DR=np.hstack((DR,dr))
        DT=np.hstack((DT,theta))
        #print len(theta)
        PSIlist[pindex]=np.sqrt(np.real(np.sum(np.exp(factor*1j*theta))/len(theta))**2+np.imag(np.sum(np.exp(factor*1j*theta))/len(theta))**2)
        pp=vor.point_region[pindex]
        polygon = [vor.vertices[i] for i in vor.regions[pp]]
        if len(neighbor_indices) >= limn and -1 not in vor.regions[pp]:
            PSI+=np.real(np.sum(np.exp(factor*1j*theta))/len(theta))
            PSII+=np.imag(np.sum(np.exp(factor*1j*theta))/len(theta))
            N+=1
    PSI=PSI/N
    PSII=PSII/N
    if visualcheck == True and ax is not None:
        voronoi_plot_2d(vor,ax=ax)
        for pindex in range(vor.npoints):
            pp=vor.point_region[pindex]
            polygon = [vor.vertices[i] for i in vor.regions[pp]]
            if -1 not in vor.regions[pp]:
                try: 
                    ax.fill(*zip(*polygon),color=mpl.cm.jet(PSIlist[pindex]))
                    #ax.annotate(str(PSIlist[pindex]),xy=(x_list[pindex], y_list[pindex]),xycoords='data')
                except: pass

    return PSI,PSII,N,PSIlist



'''
fig1=plt.figure(1,figsize=(20,10))
ax1a = fig1.add_subplot(231)
ax1b = fig1.add_subplot(234, projection='polar')
ax2a = fig1.add_subplot(232)
ax2b = fig1.add_subplot(235, projection='polar')
ax3a = fig1.add_subplot(233)
ax3b = fig1.add_subplot(236, projection='polar')

#### 1...   CHECK with real regular array ####
xy=np.meshgrid([0,2,3,5,6],[0,2,3,5,6])
x_list=np.reshape(xy[0],25)
x_list=np.append(x_list,[1,1,4,4])
y_list=np.reshape(xy[1],25)
y_list=np.append(y_list,[1,4,1,4])
xy=np.vstack((x_list, y_list)).transpose()
print "--------------------",x_list

ax1a.plot(x_list,y_list,linestyle="none",marker="o")
angular_correlation(xy,ax=ax1b)
PSI,PSII=angular_order_parameter(x_list,y_list,visualcheck=False)
print PSI,PSII
ax1a.set_xticks([])
ax1a.set_yticks([])
ax1a.set_axes("equal")

ax1b.set_xticks([0,np.pi/2.,np.pi,3*np.pi/2.])
ax1b.set_yticks([])

chain='test, $\\Psi=$%0.2f' % np.sqrt(PSI**2+PSII**2)
ax1a.annotate(chain,xy=(.5, 1),xycoords='axes fraction',xytext=(0, 20), textcoords='offset points',size=30,horizontalalignment='center',verticalalignment='bottom')


######## 2. H III H ########
data=np.loadtxt("PosH3.xls",skiprows=1)
xy=data[:,1:]
x_list=xy[:,0]
y_list=xy[:,1]

ax2a.plot(x_list,y_list,linestyle="none",marker="o")
angular_correlation(xy,ax=ax2b)
PSI,PSII=angular_order_parameter(x_list,y_list,visualcheck=False)
print PSI,PSII
ax2a.set_xticks([])
ax2a.set_yticks([])
ax2a.set_axes("equal")

ax2b.set_xticks([0,np.pi/2.,np.pi,3*np.pi/2.])
ax2b.set_yticks([])

chain='Hierar. Links, $\\Psi=$%0.2f' % np.sqrt(PSI**2+PSII**2)
ax2a.annotate(chain,xy=(.5, 1),xycoords='axes fraction',xytext=(0, 20), textcoords='offset points',size=30,horizontalalignment='center',verticalalignment='bottom')

######## Equilink ########
data=np.loadtxt("PosEQUI.xls",skiprows=1)
xy=data[:,1:]
x_list=xy[:,0]
y_list=xy[:,1]

ax3a.plot(x_list,y_list,linestyle="none",marker="o")
angular_correlation(xy,ax=ax3b)
PSI,PSII=angular_order_parameter(x_list,y_list,visualcheck=False)
print PSI,PSII
ax3a.set_xticks([])
ax3a.set_yticks([])
ax3a.set_axes("equal")

ax3b.set_xticks([0,np.pi/2.,np.pi,3*np.pi/2.])
ax3b.set_yticks([])

chain='Equilinks, $\\Psi=$%0.2f' % np.sqrt(PSI**2+PSII**2)
ax3a.annotate(chain,xy=(.5, 1),xycoords='axes fraction',xytext=(0, 20), textcoords='offset points',size=30,horizontalalignment='center',verticalalignment='bottom')


#plt.tight_layout()
fig1.savefig("order_parameter.png")
#fig1.savefig("order_parameter.eps")
'''