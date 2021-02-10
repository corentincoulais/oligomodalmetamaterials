# Save by David on 2019_07_11-12.10.11; build 2018 2017_11_07-18.21.41 127140
# This model allows for changing the thickness of hinges
# FINAL PRODUCTION: SCALE EVERYTHING BY 0.8. ALSO: Export Square hinges without HolesSquare
Mdb()

# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


import numpy as np
import os
import time

from numpy import sin, cos, tan, sqrt, pi

Model_version = 13
Subversion = 40

FindAtOn = False
AssembleHinges = True
AssembleTriangularHinges = True
steptype ='ViscoStep'	# 'Buckling' or 'StaticNL' or 'ViscoStep'
LoadingWay='DispImp' # 'Force' or 'DispImp'
if steptype=='Buckling':
	LoadY = -1.
	DispY = -1.
else:
	LoadY = -20.
	DispY = -12.5*1.5

Imperfload = 0.00 #diagonal imperfection (percentage)

ImperfectDeg =2.0#0.7 #Maximum imperfection (only if order_hinge=7)
UseStabilisation = False  # Using automatic stabilisation

SubmitJob=False
writeInput = True
ImperfGrav = False
PrecisionLevelExp = DOUBLE # SINGLE or DOUBLE
UseConstraints = False   #Applies to unrotated only!
UseConstraintsRotated = True # Applies to rotated only!
HolesSquare=False	#True if production with holes square and pocket underneath for walls. However, hinges cuts will go wrong
MergeAssembly=True
WallsMold = False
MergeIndividualParts=True 
CutHingesSquare = True #MergeIndividualParts must be true
if not MergeIndividualParts:
	CutHingesSquare =False
Stretch_fac_HingesSq = 1.05 #Lengthen square hinges for meshing (only if MergeIndividualParts)	(typically 1.04 for 1 deg deformation without wrat_hinge
	
CutwTriangularHinge = False
RotatedOrientation=True
if RotatedOrientation:
	UseConstraints = False   #Applies to unrotated only!
	
WriteSTP = False

TypeHingeTriangle = 0 #0=final geometry which can be analysed (ONLY analysis digit), 1 = geometry with connections for casting triangle hinge
MechanicalConnectionsTri = False  #To cut holes in triangular hinge for mechanical connection casting


Solvetype = 'Implicit' #'Explicit or 'Implicit'
#Time_exp=1.
MinSteps = 300
linbulk = 0.06
quadbulk = 1.2
alpha_rigid = 0.1
alpha_contact = 0.1

if TypeHingeTriangle ==0:
	ModelDim = 2	#2D Solid
	CreateMoldAssembly = False
else: 
	ModelDim = 3	#3D Solid
	Extrudedepth = 20.
	CreateMoldAssembly = True
	
#file_location =r"C:\Users\David\Documents\0_PhD\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v0"+str(Model_version)
file_location =r"F:\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_40_2deg_Visco-Agilus_improved_norelax19mm"
file_name = file_location+'\\Model_v0'+str(Model_version)+'_0'+str(Subversion)+'.cae'
os.chdir( file_location)


if steptype =='Buckling':
	jobtitlestart = 'Job-Buck_v0'+str(Model_version)+'_0'+str(Subversion)	
elif steptype =='StaticNL':
	jobtitlestart = 'Job-NLCompression_v0'+str(Model_version)+'_0'+str(Subversion)	
elif steptype == 'ViscoStep':
	jobtitlestart = 'Job-NLViscoComp_v0'+str(Model_version)+'_0'+str(Subversion)	
def maketuple(dic):
	a=[]
	for ind in range(len(dic)):
		a=a+[dic[dic.keys()[ind]]]
	a=tuple(a)
	return a


BCTop = 'support' # 'rigid' or 'support'
BCBottom = 'support' # 'rigid' or 'support'

if BCTop =='rigid' or BCBottom =='rigid':
	UseContact = True
else:
	UseContact = False
FricCoef = 0.05
RigidMass = 0.1


amptype ='Smooth' #'Straight'  (linear) or 'Smooth'



#########################
#Materials
#########################

E_Rigid = 1000.
nu_Rigid = 0.3
E_SqH=  1.
E_TriH_vec = 6*[1.]+6*[2.]+6*[0.5]		+[1.,2.,0.5]
#E_TriH_vec = 6*[10.]+6*[5.]+6*[2.]


nu_H = 0.47
#DensMat = 1e-4  #vector now for stability

MatPot = 'ArrudaBoyceVisco' #'NeoHooke' Or 'Linear' or 'ArrudaBoyce' or 'ArrudaBoyceVisco' 
lambda_SqH = E_SqH*nu_H/((1+nu_H)*(1-2*nu_H))
mu_SqH=E_SqH/(2*(1+nu_H))
C10_SqH = mu_SqH/2
D1_SqH = lambda_SqH/2

lam_AB = 2 # #stiffenening parameter arruda-boyce

Pronytermsvisco = len(E_TriH_vec)*[((0.4472617812538507, 0.4472617812538507, 0.04742040029399319),(0.2557510036456948, 0.2557510036456948, 0.972313952967699),(0.10196933865360647, 0.10196933865360647, 17.783982821444095))]

#Damping (explicit only)
alpha_damp = [0.05,0.05,0.05]#[0.05,0.05,0.05]
beta_damp = [0.05,0.05,0.05]

ind_ItMat = 0
N_ItMat = len(E_TriH_vec)



#########################
# Relaxation and time specs
#########################

t_loading_vec = 3*[1.5e-1,1.5e0,1.5e1,1.5e2,1.5e3,1.5e4]		+3*[1.5e5]
t_loading=t_loading_vec[0]

density_vec = 3*[1e-3,1e-2,1e-1,1e0,1e1,1e2]		+3*[1.e3]
DensMat = density_vec[0]

Dorelaxation = False



t_relaxstart  =0.
t_relaxend = 1e3
dt_relaxstart=1e-2
facdtrelax = 1.1

facrangerelax=np.array(range(0,1000))
dt_relaxrange = dt_relaxstart*facdtrelax**facrangerelax
t_relaxlong = t_relaxstart+np.cumsum(dt_relaxrange)
t_relax=t_relaxlong[:np.where(t_relaxlong>t_relaxend)[0][0]]
lentrelax = len(t_relax)
t_relaxtuple=[]
for indrelax in range(lentrelax):
	t_relaxtuple+=[(1e-4*np.round(1e4*t_relax[indrelax]),)]
t_relaxtuple=tuple(t_relaxtuple)




#########################
#Geometry
#########################

#Geometry main  unit cell
a_cell = 15. #Unit cell size
d_cell =  np.sqrt(2)*a_cell #Unit cell diagonal
phi_hinge_deg = 25. #25. for everything before
phi_hinge =  np.deg2rad(phi_hinge_deg)
phi_hinge_deg_tri = 25. #25. for everything before
phi_hinge_tri =  np.deg2rad(phi_hinge_deg_tri)
wrat_hinge = 0.95 #ratio  #0.95 in final sample
wrat_hinge_tri = 0.99  #0.99 in final sample
w_hinge = 4./wrat_hinge  #4./wrat_hinge in final sample
w_hinge_tri = 4./wrat_hinge_tri #4./wrat_hinge_tri in final sample
r_fillet=0.17#.2/wrat_hinge_tri in final sample
r_fillet_tri = 0.17
s_tri = 0.3





t_edge_break  = 0.05
t_max_break = 0.5


dh = d_cell/2
wh = w_hinge/2 #Outer hinge on solid
hh = wh/np.tan(phi_hinge)
bh = dh-hh
lh = wh/np.sin(phi_hinge)

wh2 = wh*wrat_hinge		#Actual hinge
lh2 = wh2/np.sin(phi_hinge)
pv2 = wh2/tan(phi_hinge)*2 


dh3 = dh-r_fillet/sin(phi_hinge)+r_fillet#distance to fillet edges

wh4 = wh2+2*t_edge_break		#Hinge after breaking, maybe edit to thickness?
lh4 = wh4/np.sin(phi_hinge)
wbmax = wh4+t_max_break



wh5 = w_hinge_tri/2 #Outer hinge on solid
hh5 = wh5/np.tan(phi_hinge_tri)
bh5 = dh-hh5
lh5 = wh5/np.sin(phi_hinge_tri)

wh6 = wh5*wrat_hinge_tri		#Actual hinge TRIANGULAR
lh6 = wh6/np.sin(phi_hinge_tri)
pv6 = wh6/tan(phi_hinge_tri)*2


#centre[0] top corner [1-3], right corner[4-6], bottom corner[7-9], right corner[10-12], 4radius centres [13-16]
xv = [0, 	0,wh,-wh, 		dh,bh,bh,		0,-wh,wh,		-dh,-bh,-bh,		dh,dh,-dh,-dh]
yv = [0, 	dh,bh,bh,		0,-wh,wh,	 	-dh,-bh,-bh,    0,wh,-wh,			dh, -dh, -dh, dh]

rv1 = np.sqrt((xv[13]-xv[6])**2+(yv[13]-yv[6])**2)
sv1 = d_cell*np.sqrt(2)
wv1 = sv1-2.*rv1
rv1 = 3.
depthv1=4.
tv1=5.

#Geometry secondary triangular unit cell
theta = np.pi/4-phi_hinge
theta_tri = np.pi/4-phi_hinge_tri
xv2=xv[:]
yv2=yv[:]
xv2=xv2+[lh5*sin(theta_tri),lh5*cos(theta_tri), -lh5*sin(theta_tri),-lh5*cos(theta_tri)] #[17-20] (edges triangular hinge holders)
yv2=yv2+[lh5*cos(theta_tri),lh5*sin(theta_tri),-lh5*cos(theta_tri),-lh5*sin(theta_tri) ] #[17-20](edges triangular hinge holders)
"""
xv2=xv2+[lh*sin(theta),lh*cos(theta), -lh*sin(theta),-lh*cos(theta)] #[17-20]
yv2=yv2+[lh*cos(theta),lh*sin(theta),-lh*cos(theta),-lh*sin(theta) ] #[17-20]
"""

xv2 = xv2+[xv2[17]+s_tri/sqrt(2), xv2[18]+s_tri/sqrt(2), xv2[19]-s_tri/sqrt(2), xv2[20]-s_tri/sqrt(2)] #[21-23]
yv2 = yv2+[yv2[17]+s_tri/sqrt(2), yv2[18]+s_tri/sqrt(2), yv2[19]-s_tri/sqrt(2), yv2[20]-s_tri/sqrt(2)] #[21-23]

xv2=xv2+[lh2*sin(theta),lh2*cos(theta), -lh2*sin(theta),-lh2*cos(theta)] #[25-28] (central hinge)
yv2=yv2+[lh2*cos(theta),lh2*sin(theta),-lh2*cos(theta),-lh2*sin(theta) ] #[25-28](central hinge)

# Geometry attachments breaking [29-34]

xv2=xv2+[lh4*sin(theta),lh4*cos(theta), -lh4*sin(theta),-lh4*cos(theta)] #[29-32] (corners break area)
yv2=yv2+[lh4*cos(theta),lh4*sin(theta),-lh4*cos(theta),-lh4*sin(theta) ] #[29-32](corners break area)
xv2=xv2+[wbmax/sqrt(2),-wbmax/sqrt(2)]#[33-34] (centres break area)
yv2=yv2+[-wbmax/sqrt(2),wbmax/sqrt(2)]#[33-34] (centres break area)


xv2=xv2+[lh6*sin(theta_tri),lh6*cos(theta_tri), -lh6*sin(theta_tri),-lh6*cos(theta_tri)] #[35-38] (central TRIANGULAR hinge)
yv2=yv2+[lh6*cos(theta_tri),lh6*sin(theta_tri),-lh6*cos(theta_tri),-lh6*sin(theta_tri) ] #[35-38](central TRIANGULAR hinge)



xv=xv2[:]
yv=yv2[:]


#######################
# Assembly
#######################

#Index 0: regular square, index 1: double triangle with hinge, index 2: double triangle witout hing but with walls? 10: NOT PRESENT

if RotatedOrientation:
	Nx=12
	Ny=12
	
	if Nx==9 and Ny==9:
		Assembly_mat=[
		[-1,-1,-1,-1,0,-1,-1,-1,-1],
		[-1,-1,-1,0,0,0,-1,-1,-1],
		[-1,-1,0,0,0,0,0,-1,-1],
		[-1,0,0,0,0,0,0,0,-1],
		[0,0,0,0,0,0,0,0,0],
		[-1,0,0,0,0,0,0,0,-1],
		[-1,-1,0,0,0,0,0,-1,-1],
		[-1,-1,-1,0,0,0,-1,-1,-1],
		[-1,-1,-1,-1,0,-1,-1,-1,-1]]
		
		Nxrot = Nx+1
		Nyrot  = Ny+1
		Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	
		
	elif Nx==10 and Ny==10:
		Assembly_mat=[
		[-1,-1,-1,-1,0,0,-1,-1,-1,-1],
		[-1,-1,-1,0,0,0,0,-1,-1,-1],
		[-1,-1,0,0,0,0,0,0,-1,-1],
		[-1,0,0,0,0,0,0,0,0,-1],
		[0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0],
		[-1,0,0,0,0,0,0,0,0,-1],
		[-1,-1,0,0,0,0,0,0,-1,-1],
		[-1,-1,-1,0,0,0,0,-1,-1,-1],
		[-1,-1,-1,-1,0,0,-1,-1,-1,-1]]
		
		Nxrot = Nx+1
		Nyrot  = Ny+1
		NBCrot = int((Nxrot-1)/2)
		DeltaBCrot = NBCrot
		Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	

		
	elif Nx==12 and Ny==12:
		orderhinge=8 #2 for model v08_01, 3 for v10_01
		if orderhinge==0:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		0,0,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  0,0,0,0,			-1,-1,-1,-1],
			[-1,-1,-1,			0,0,0,0,0,0,			-1,-1,-1],
			[-1,-1,			  0,0,0,0,0,0,0,0,			-1,-1],
			[-1, 		    0,0,0,0,0,0,0,0,0,0,			-1],
			[			  0,0,0,0,0,0,0,0,0,0,0,0],
			[			  0,0,0,0,0,0,0,0,0,0,0,0],
			[-1,            0,0,0,0,0,0,0,0,0,0,      	-1],
			[-1,-1,           0,0,0,0,0,0,0,0,     		-1,-1],
			[-1,-1,-1,          0,0,0,0,0,0,        	-1,-1,-1],
			[-1,-1,-1,-1,         0,0,0,0,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        0,0, 				-1,-1,-1,-1,-1]]
				
			Nxrot = Nx+1
			Nyrot  = Ny+1
			NBCrot = int((Nxrot-1)/2)
			DeltaBCrot = NBCrot
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	

		elif orderhinge==1:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		0,1,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  0,0,0,0,			-1,-1,-1,-1],
			[-1,-1,-1,			0,1,0,1,0,1,			-1,-1,-1],
			[-1,-1,			  0,0,0,0,0,0,0,0,			-1,-1],
			[-1, 		    0,1,0,1,0,1,0,1,0,1,			-1],
			[			  0,0,0,0,0,0,0,0,0,0,0,0		],
			[			  1,0,1,0,1,0,1,0,1,0,1,0		],
			[-1,            0,0,0,0,0,0,0,0,0,0,      	-1],
			[-1,-1,           1,0,1,0,1,0,1,0,     		-1,-1],
			[-1,-1,-1,          0,0,0,0,0,0,        	-1,-1,-1],
			[-1,-1,-1,-1,         1,0,1,0,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        0,0, 				-1,-1,-1,-1,-1]]
				
			Nxrot = Nx+1
			Nyrot  = Ny+1
			NBCrot = int((Nxrot-1)/2)
			DeltaBCrot = NBCrot
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	

		elif orderhinge==2:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		0,0,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  1,0,1,0,			-1,-1,-1,-1],
			[-1,-1,-1,			0,0,0,0,0,0,			-1,-1,-1],
			[-1,-1,			  1,0,1,0,1,0,1,0,			-1,-1],
			[-1, 		    0,0,0,0,0,0,0,0,0,0,			-1],
			[			  1,0,1,0,1,0,1,0,1,0,1,0		],
			[			  0,0,0,0,0,0,0,0,0,0,0,0		],
			[-1,            0,1,0,1,0,1,0,1,0,1,      	-1],
			[-1,-1,           0,0,0,0,0,0,0,0,     		-1,-1],
			[-1,-1,-1,          0,1,0,1,0,1,        	-1,-1,-1],
			[-1,-1,-1,-1,         0,0,0,0,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        0,1, 				-1,-1,-1,-1,-1]]	
				
			Nxrot = Nx+1
			Nyrot  = Ny+1
			NBCrot = int((Nxrot-1)/2)
			DeltaBCrot = NBCrot
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	


		elif orderhinge==3:
			Assembly_mat=[			
			[-1,-1,-1,-1,-1,		0,0,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  1,0,1,0,			-1,-1,-1,-1],
			[-1,-1,-1,			0,0,0,0,0,0,			-1,-1,-1],
			[-1,-1,			  1,0,1,0,1,0,1,0,			-1,-1],
			[-1, 		    0,0,0,0,0,0,0,0,0,0,			-1],
			[			  1,0,1,0,1,0,1,0,1,0,1,0		],
			[			  0,0,0,0,0,0,0,0,0,0,0,		-1		],
			[-1,            0,1,0,1,0,1,0,1,0,			-1,      	-1],
			[-1,-1,           0,0,0,0,0,0,0,			-1,     		-1,-1],
			[-1,-1,-1,          0,1,0,1,0,				-1,        	-1,-1,-1],
			[-1,-1,-1,-1,         0,0,0,				-1,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        0,					-1, 				-1,-1,-1,-1,-1]]	
			Nxrot = 13#Nx+1
			Nyrot  = 13#Ny+1
			NBCrot = 6#int((Nxrot-1)/2)
			DeltaBCrot = NBCrot
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	
			

		elif orderhinge==4:
			Assembly_mat=[			
			[-1,-1,-1,-1,-1,		0,0,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  1,0,1,0,			-1,-1,-1,-1],
			[-1,-1,-1,			0,0,0,0,0,0,			-1,-1,-1],
			[-1,-1,			  1,0,1,0,1,0,1,0,			-1,-1],
			[-1, 		    0,0,0,0,0,0,0,0,0,-1,			-1],
			[			  1,0,1,0,1,0,1,0,1,		-1,-1,		-1		],
			[			  0,0,0,0,0,0,0,0,		-1,-1,		-1,		-1		],
			[-1,            0,1,0,1,0,1,		-1,-1,		-1,			-1,      	-1],
			[-1,-1,           0,0,0,0,		-1,-1,		-1,			-1,     		-1,-1],
			[-1,-1,-1,          0,1,		-1,-1,		-1,				-1,        	-1,-1,-1],
			[-1,-1,-1,-1,         		-1,-1,		-1,				-1,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        	-1,					-1, 				-1,-1,-1,-1,-1]]				
							
			Nxrot = 13#Nx+1
			Nyrot  = 13#Ny+1
			NBCrot = 4#int((Nxrot-1)/2)
			DeltaBCrot = 6
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	
			
		elif orderhinge==7:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		3,3,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  3,0,1,3,			-1,-1,-1,-1],
			[-1,-1,-1,			3,0,0,0,0,3,			-1,-1,-1],
			[-1,-1,			  3,0,1,0,1,0,1,3,			-1,-1],
			[-1, 		    3,0,0,0,0,0,0,0,0,3,			-1],
			[			  3,0,1,0,1,0,1,0,1,0,1,3		],
			[			  3,0,0,0,0,0,0,0,0,0,0,3		],
			[-1,            3,1,0,1,0,1,0,1,0,3,      	-1],
			[-1,-1,           3,0,0,0,0,0,0,3,     		-1,-1],
			[-1,-1,-1,          3,1,0,1,0,3,        	-1,-1,-1],
			[-1,-1,-1,-1,         3,0,0,3,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        3,3, 				-1,-1,-1,-1,-1]]	
				
			Nxrot = Nx+1
			Nyrot  = Ny+1
			NBCrot = int((Nxrot-1)/2)
			DeltaBCrot = NBCrot
			#Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated		
			a=np.float(ImperfectDeg)/6 #Deg change per time
			Rotate_mat=	[
			[0,0,0,0,0,			    6,-5,			0,0,0,0,0],
			[0,0,0,0,		 	  6,-5,4,-3,			0,0,0,0],
			[0,0,0,				6,-5,4,-3,2,-1,			0,0,0],
			[0,0,			  6,-5,4,-3,2,-1,0,1,			0,0],
			[0, 		    6,-5,4,-3,2,-1,0,1,-2,3,			0],
			[			  6,-5,4,-3,2,-1,0,1,-2,3,-4,5		],
			[			  -5,4,-3,2,-1,0,1,-2,3,-4,5,-6		],
			[0,              -3,2,-1,0,1,-2,3,-4,5,-6,      	0],
			[0,0,               -1,0,1,-2,3,-4,5,-6,     		0,0],
			[0,0,0,                1,-2,3,-4,5,-6,        	0,0,0],
			[0,0,0,0,                 3,-4,5,-6,          	0,0,0,0],
			[0,0,0,0,0,                  5,-6, 				0,0,0,0,0]]
			
			Rotate_mat_mainonly=np.array([
			[0.,0.,0.,0.,0.,		-90,180,			0.,0.,0.,0.,0.],
			[0.,0.,0.,0.,		  -90,0,0,180,			0.,0.,0.,0.],
			[0.,0.,0.,			-90,0,0,0,0,180,			0.,0.,0.],
			[0.,0.,			  -90,0,0,0,0,0,0,180,			0.,0.],
			[0., 		    -90,0,0,0,0,0,0,0,0,180,			0.],
			[			  -90,0,0,0,0,0,0,0,0,0,0,180],
			[			  0,0,0,0,0,0,0,0,0,0,0,90],
			[0.,            0,0,0,0,0,0,0,0,0,90,      	0.],
			[0.,0.,           0,0,0,0,0,0,0,90,     		0.,0.],
			[0.,0.,0.,          0,0,0,0,0,90,        	0.,0.,0.],
			[0.,0.,0.,0.,         0,0,0,90,          	0.,0.,0.,0.],
			[0.,0.,0.,0.,0.,        0,90, 				0.,0.,0.,0.,0.]],dtype=float)		
			"""
			Rotate_mat=[			
			[0.,0.,0.,0.,0.,		6.*a,-5.*a,			0.,0.,0.,0.,0.],
			[0.,0.,0.,0.,		  6.*a,-5.*a,4.*a,-3.*a,			0.,0.,0.,0.],
			[0.,0.,0.,			6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,			0.,0.,0.],
			[0.,0.,			  6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,			0.,0.],
			[0., 		    6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,			0.],
			[			  6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a		],
			[			  -5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a,-6.*a		],
			[0.,            -3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a,-6.*a,      	0.],
			[0.,0.,           -1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a,-6.*a,     		0.,0.],
			[0.,0.,0.,          -1.*a,2.*a,-3.*a,4.*a,5.*a,-6.*a,        	0.,0.,0.],
			[0.,0.,0.,0.,         3.*a,-4.*a,5.*a,-6.*a,          	0.,0.,0.,0.],
			[0.,0.,0.,0.,0.,        5.*a,-6.*a, 				0.,0.,0.,0.,0.]]	
			"""
			Rotate_mat=np.array(Rotate_mat,dtype=float)
			Rotate_mat=a*Rotate_mat
			
		elif orderhinge==8:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		3,3,			-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  3,0,0,3,			-1,-1,-1,-1],
			[-1,-1,-1,			3,1,0,1,0,3,			-1,-1,-1],
			[-1,-1,			  3,0,0,0,0,0,0,3,			-1,-1],
			[-1, 		    3,1,0,1,0,1,0,1,0,3,			-1],
			[			  3,0,0,0,0,0,0,0,0,0,0,3		],
			[			  3,0,1,0,1,0,1,0,1,0,1,3		],
			[-1,            3,0,0,0,0,0,0,0,0,3,      	-1],
			[-1,-1,           3,0,1,0,1,0,1,3,     		-1,-1],
			[-1,-1,-1,          3,0,0,0,0,3,        	-1,-1,-1],
			[-1,-1,-1,-1,         3,0,1,3,          	-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,        3,3, 				-1,-1,-1,-1,-1]]

			
			for ind1 in range(len(Assembly_mat)):
				for ind2 in range(len(Assembly_mat[0])):
					if Assembly_mat[ind1][ind2]==3: Assembly_mat[ind1][ind2]=0
			
			#Assembly_mat = np.transpose(Assembly_mat)
			bla=1
			Nxrot = Nx+1
			Nyrot  = Ny+1
			NBCrot = int((Nxrot-1)/2)
			DeltaBCrot = NBCrot
			#Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated		
			a=np.float(ImperfectDeg)/6. #Deg change per time
			Rotate_mat=	[
			[0,0,0,0,0,			    6,-5,			0,0,0,0,0],
			[0,0,0,0,		 	  6,-5,4,-3,			0,0,0,0],
			[0,0,0,				6,-5,4,-3,2,-1,			0,0,0],
			[0,0,			  6,-5,4,-3,2,-1,0,1,			0,0],
			[0, 		    6,-5,4,-3,2,-1,0,1,-2,3,			0],
			[			  6,-5,4,-3,2,-1,0,1,-2,3,-4,5		],
			[			  -5,4,-3,2,-1,0,1,-2,3,-4,5,-6		],
			[0,              -3,2,-1,0,1,-2,3,-4,5,-6,      	0],
			[0,0,               -1,0,1,-2,3,-4,5,-6,     		0,0],
			[0,0,0,                1,-2,3,-4,5,-6,        	0,0,0],
			[0,0,0,0,                 3,-4,5,-6,          	0,0,0,0],
			[0,0,0,0,0,                  5,-6, 				0,0,0,0,0]]
			
			Rotate_mat = np.array(Rotate_mat,dtype=float)[:,::-1]#np.transpose(Rotate_mat)
			
			
			Rotate_mat_mainonly=np.array([
			[0.,0.,0.,0.,0.,		-90,180,			0.,0.,0.,0.,0.],
			[0.,0.,0.,0.,		  -90,0,0,180,			0.,0.,0.,0.],
			[0.,0.,0.,			-90,0,0,0,0,180,			0.,0.,0.],
			[0.,0.,			  -90,0,0,0,0,0,0,180,			0.,0.],
			[0., 		    -90,0,0,0,0,0,0,0,0,180,			0.],
			[			  -90,0,0,0,0,0,0,0,0,0,0,180],
			[			  0,0,0,0,0,0,0,0,0,0,0,90],
			[0.,            0,0,0,0,0,0,0,0,0,90,      	0.],
			[0.,0.,           0,0,0,0,0,0,0,90,     		0.,0.],
			[0.,0.,0.,          0,0,0,0,0,90,        	0.,0.,0.],
			[0.,0.,0.,0.,         0,0,0,90,          	0.,0.,0.,0.],
			[0.,0.,0.,0.,0.,        0,90, 				0.,0.,0.,0.,0.]],dtype=float)		
			

			for ind1 in range(len(Assembly_mat)):
				for ind2 in range(len(Assembly_mat[0])):
					if Assembly_mat[ind1][ind2]==1: Rotate_mat_mainonly[ind1][ind2]=90			
			"""
			Rotate_mat=[			
			[0.,0.,0.,0.,0.,		6.*a,-5.*a,			0.,0.,0.,0.,0.],
			[0.,0.,0.,0.,		  6.*a,-5.*a,4.*a,-3.*a,			0.,0.,0.,0.],
			[0.,0.,0.,			6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,			0.,0.,0.],
			[0.,0.,			  6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,			0.,0.],
			[0., 		    6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,			0.],
			[			  6.*a,-5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a		],
			[			  -5.*a,4.*a,-3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a,-6.*a		],
			[0.,            -3.*a,2.*a,-1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a,-6.*a,      	0.],
			[0.,0.,           -1.*a,0.*a,1.*a,-2.*a,3.*a,-4.*a,5.*a,-6.*a,     		0.,0.],
			[0.,0.,0.,          -1.*a,2.*a,-3.*a,4.*a,5.*a,-6.*a,        	0.,0.,0.],
			[0.,0.,0.,0.,         3.*a,-4.*a,5.*a,-6.*a,          	0.,0.,0.,0.],
			[0.,0.,0.,0.,0.,        5.*a,-6.*a, 				0.,0.,0.,0.,0.]]	
			"""
			Rotate_mat=np.array(Rotate_mat,dtype=float)
			Rotate_mat=a*Rotate_mat
			
						

	elif Nx==13 and Ny==13:
		orderhinge=5 #2 for model v08_01, 3 for v10_01

		if orderhinge==5:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		0,0,			-1,-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  1,0,1,0,			-1,-1,-1,-1,-1],
			[-1,-1,-1,			0,0,0,0,0,0,			-1,-1,-1,-1],
			[-1,-1,			  1,0,1,0,1,0,1,0,			-1,-1,-1],
			[-1, 		    0,0,0,0,0,0,0,0,0,0,			-1,-1],
			[			  1,0,1,0,1,0,1,0,1,0,1,0		,-1],
			[			  0,0,0,0,0,0,0,0,0,0,0,0,0		],
			[-1,            0,1,0,1,0,1,0,1,0,1,0,1      	],
			[-1,-1,           0,0,0,0,0,0,0,0,0,0,     		-1],
			[-1,-1,-1,          0,1,0,1,0,1,0,1,        	-1,-1],
			[-1,-1,-1,-1,         0,0,0,0,0,0,          	-1,-1,-1],
			[-1,-1,-1,-1,-1,        0,1,0,1, 				-1,-1,-1,-1],	
			[-1,-1,-1,-1,-1,-1,       0,0,      			-1,-1,-1,-1,-1]]
				
			Nxrot = Nx+1
			Nyrot  = 13#Ny+1
			NBCrot = 7#int((Nxrot-1)/2)
			DeltaBCrot = 6	
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated	
	elif Nx==14 and Ny==14:
		orderhinge=6 #2 for model v08_01, 3 for v10_01

		if orderhinge==6:
			Assembly_mat=[
			[-1,-1,-1,-1,-1,		0,0,			-1,-1,-1,-1,-1,-1,-1],
			[-1,-1,-1,-1,		  1,0,1,0,			-1,-1,-1,-1,-1,-1],
			[-1,-1,-1,			0,0,0,0,0,0,			-1,-1,-1,-1,-1],
			[-1,-1,			  1,0,1,0,1,0,1,0,			-1,-1,-1,-1],
			[-1, 		    0,0,0,0,0,0,0,0,0,0,			-1,-1,-1],
			[			  1,0,1,0,1,0,1,0,1,0,1,0,       -1,-1],
			[			  0,0,0,0,0,0,0,0,0,0,0,0,0,		-1],
			[-1,            0,1,0,1,0,1,0,1,0,1,0,1,0      	],
			[-1,-1,           0,0,0,0,0,0,0,0,0,0,0,0     	],
			[-1,-1,-1,          0,1,0,1,0,1,0,1,0,1,        	-1],
			[-1,-1,-1,-1,         0,0,0,0,0,0,0,0,          	-1,-1],
			[-1,-1,-1,-1,-1,        0,1,0,1,0,1, 				-1,-1,-1],	
			[-1,-1,-1,-1,-1,-1,       0,0,0,0,      			-1,-1,-1,-1],
			[-1,-1,-1,-1,-1,-1,-1,      0,1,      			-1,-1,-1,-1,-1]]
				
			Nxrot = Nx+1
			Nyrot  = 13#Ny+1
			NBCrot = 8#int((Nxrot-1)/2)
			DeltaBCrot = 6			
			Rotate_mat=0.*np.ones((Nx,Ny))   #Determines how samples are rotated				
	Assembly_mat=np.array(Assembly_mat)
	
else:
	Nx = 7
	Ny = 7

	Assembly_mat = np.zeros((Nx,Ny))

	for indX in [1,3,5]:
		for indY in [1,3,5]:
			Assembly_mat[indX,indY] = 1

	Assembly_mat = Assembly_mat.astype(int)

h_rigid = 1. #Distance above sample initially
h_rigid_RP = 0. #distance rotation point rigid above rigid

### MESHSIZE
MeshSize = d_cell/8.
SeedEdge = True
SeedCountSide =2
SeedCountStraight =3
SeedCountHinge =1

GravLoad = 0.
FixX = False
UseLoadY = True

K_Contact =5.#Contact stiffness










#########################
#Model
#########################

Mdb()

# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


if FindAtOn:
	session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
else:
	session.journalOptions.setValues(replayGeometry=INDEX, recoverGeometry=INDEX)

#########################
# Create Materials and Sections
#########################

Materials = ['Rigid','HingeSq','HingeTri']
SectionsShell = ['Rigid_Shell','HingeSq_Shell','HingeTri_Shell']
SectionsSol = ['Rigid_Solid','HingeSq_Solid','HingeTri_Solid']

lenSecSol = len(SectionsSol)
SectionsAll = SectionsSol+SectionsShell
if ModelDim==3:
	dsec = 0#lenSecSol
else:
	dsec = 0
	

Mat  = Materials[0]
Sec = SectionsShell[0]
SecSol = SectionsSol[0]

mdb.models['Model-1'].Material(name=Mat)
mdb.models['Model-1'].materials[Mat].Density(table=((DensMat, ), ))
mdb.models['Model-1'].materials[Mat].Elastic(table=((E_Rigid, nu_Rigid), ))
if Solvetype =='Explicit':
	mdb.models['Model-1'].materials[Mat].Damping(alpha=alpha_damp[0],beta=beta_damp[0])

mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material=Mat, name=Sec, 
    nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    preIntegrate=OFF, temperature=GRADIENT, thickness=1.0, thicknessField='', 
    thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
mdb.models['Model-1'].HomogeneousSolidSection(material=Mat, name=
    SecSol, thickness=None)	

Mat  = Materials[1]
Sec = SectionsShell[1]
SecSol = SectionsSol[1]

mdb.models['Model-1'].Material(name=Mat)
mdb.models['Model-1'].materials[Mat].Density(table=((DensMat, ), ))
if MatPot =='Linear':
	mdb.models['Model-1'].materials[Mat].Elastic(table=((E_SqH, nu_H), ))
elif MatPot == 'NeoHooke':	
	mdb.models['Model-1'].materials[Mat].Hyperelastic(materialType=
		ISOTROPIC, table=((C10_SqH, D1_SqH), ), testData=OFF, type=NEO_HOOKE, 
		volumetricResponse=VOLUMETRIC_DATA)
elif MatPot =='ArrudaBoyce' :
	mdb.models['Model-1'].materials[Mat].Hyperelastic(materialType=
		ISOTROPIC, table=((2*C10_SqH, lam_AB, D1_SqH), ), testData=OFF, type=ARRUDA_BOYCE, 
		volumetricResponse=VOLUMETRIC_DATA)
elif MatPot =='ArrudaBoyceVisco':
	mdb.models['Model-1'].materials[Mat].Hyperelastic(materialType=
		ISOTROPIC, table=((2*C10_SqH, lam_AB, D1_SqH), ), testData=OFF, type=ARRUDA_BOYCE, 
		volumetricResponse=VOLUMETRIC_DATA)	
	mdb.models['Model-1'].materials[Mat].Viscoelastic(domain=TIME, table=Pronytermsvisco[0], time=PRONY)
		
		
if Solvetype =='Explicit':
	mdb.models['Model-1'].materials[Mat].Damping(alpha=alpha_damp[1],beta=beta_damp[1])

mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material=Mat, name=Sec, 
    nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    preIntegrate=OFF, temperature=GRADIENT, thickness=1.0, thicknessField='', 
    thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
mdb.models['Model-1'].HomogeneousSolidSection(material=Mat, name=
    SecSol, thickness=None)		
	
Mat  = Materials[2]
Sec = SectionsShell[2]
SecSol = SectionsSol[2]

E_TriH = E_TriH_vec[0]
lambda_TriH = E_TriH*nu_H/((1+nu_H)*(1-2*nu_H))
mu_TriH=E_TriH/(2*(1+nu_H))
C10_TriH = mu_TriH/2
D1_TriH = lambda_TriH/2

mdb.models['Model-1'].Material(name=Mat)
mdb.models['Model-1'].materials[Mat].Density(table=((DensMat, ), ))
if MatPot =='Linear':
	mdb.models['Model-1'].materials[Mat].Elastic(table=((E_TriH, nu_H), ))
elif MatPot == 'NeoHooke':	
	mdb.models['Model-1'].materials[Mat].Hyperelastic(materialType=
		ISOTROPIC, table=((C10_TriH, D1_TriH), ), testData=OFF, type=NEO_HOOKE, 
		volumetricResponse=VOLUMETRIC_DATA)
elif MatPot =='ArrudaBoyce' or MatPot =='ArrudaBoyceVisco':
	mdb.models['Model-1'].materials[Mat].Hyperelastic(materialType=
		ISOTROPIC, table=((2*C10_TriH, lam_AB, D1_TriH), ), testData=OFF, type=ARRUDA_BOYCE, 
		volumetricResponse=VOLUMETRIC_DATA)	
"""
elif MatPot =='ArrudaBoyceVisco':
	mdb.models['Model-1'].materials[Mat].Hyperelastic(materialType=
		ISOTROPIC, table=((2*C10_TriH, lam_AB, D1_TriH), ), testData=OFF, type=ARRUDA_BOYCE, 
		volumetricResponse=VOLUMETRIC_DATA)	
	mdb.models['Model-1'].materials[Mat].Viscoelastic(domain=TIME, table=(Pronytermsvisco[0], ), time=PRONY)
"""		

		
if Solvetype =='Explicit':
	mdb.models['Model-1'].materials[Mat].Damping(alpha=alpha_damp[2],beta=beta_damp[2])

mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material=Mat, name=Sec, 
    nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    preIntegrate=OFF, temperature=GRADIENT, thickness=1.0, thicknessField='', 
    thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
mdb.models['Model-1'].HomogeneousSolidSection(material=Mat, name=
    SecSol, thickness=None)		
	


#########################
#Part, Unit cell
#########################
Partname = 'UnitCellMain'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
#create corners

for ind in range(4):
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[1+3*ind], yv[1+3*ind]), point2=(
		xv[2+3*ind], yv[2+3*ind]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[1+3*ind], yv[1+3*ind]), point2=(
		xv[3+3*ind], yv[3+3*ind]))
	
	
	
	
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(xv[13], yv[13])
	, direction=COUNTERCLOCKWISE, point1=(xv[2], yv[2]), point2=(xv[6], yv[6]))
	
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(xv[14], yv[14])
	, direction=COUNTERCLOCKWISE, point1=(xv[5], yv[5]), point2=(xv[9], yv[9]))
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(xv[15], yv[15])
	, direction=COUNTERCLOCKWISE, point1=(xv[8], yv[8]), point2=(xv[12], yv[12]))
	
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(xv[16], yv[16])
	, direction=COUNTERCLOCKWISE, point1=(xv[11], yv[11]), point2=(xv[3], yv[3]))
	
idP=1
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
    xv[idP]*0.999, yv[idP]*0.999), nearPoint2=(xv[idP]*0.999, 
    yv[idP]*0.999), radius=r_fillet)
	
idP=4
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[4], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[5], nearPoint1=(
    xv[idP]*0.999, yv[idP]*0.999), nearPoint2=(xv[idP]*0.999, 
    yv[idP]*0.999), radius=r_fillet)	
		
idP=7
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[6], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[7], nearPoint1=(
    xv[idP]*0.999, yv[idP]*0.999), nearPoint2=(xv[idP]*0.999, 
    yv[idP]*0.999), radius=r_fillet)	
	
idP=10
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[8], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[9], nearPoint1=(
    xv[idP]*0.999, yv[idP]*0.999), nearPoint2=(xv[idP]*0.999, 
    yv[idP]*0.999), radius=r_fillet)	
	
	
if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].parts['UnitCellMain'].PartitionEdgeByParam(edges=
    mdb.models['Model-1'].parts['UnitCellMain'].edges[0:1]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[4:5]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[8:9]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[12:13], parameter=0.5)
	
mdb.models['Model-1'].parts['UnitCellMain'].Set(name='SquareTop', vertices=
    mdb.models['Model-1'].parts['UnitCellMain'].vertices[1:2])
mdb.models['Model-1'].parts['UnitCellMain'].Set(name='SquareRight', vertices=
    mdb.models['Model-1'].parts['UnitCellMain'].vertices[16:17])
mdb.models['Model-1'].parts['UnitCellMain'].Set(name='SquareBottom', vertices=
    mdb.models['Model-1'].parts['UnitCellMain'].vertices[11:12])
mdb.models['Model-1'].parts['UnitCellMain'].Set(name='SquareLeft', vertices=
    mdb.models['Model-1'].parts['UnitCellMain'].vertices[6:7])	


mdb.models['Model-1'].parts['UnitCellMain'].Set(edges=
    mdb.models['Model-1'].parts['UnitCellMain'].edges[3:4]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[8:9]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[13:14]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[18:19], name='SideEdges')
mdb.models['Model-1'].parts['UnitCellMain'].Set(edges=
    mdb.models['Model-1'].parts['UnitCellMain'].edges[2:3]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[4:5]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[7:8]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[9:10]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[12:13]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[14:15]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[17:18]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[19:20], name=
    'CornerStraightEdges')

mdb.models['Model-1'].parts['UnitCellMain'].Set(edges=
    mdb.models['Model-1'].parts['UnitCellMain'].edges[0:2]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[5:7]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[10:12]+\
    mdb.models['Model-1'].parts['UnitCellMain'].edges[15:17], name=
    'HingeEdges')


	
#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[0+dsec], thicknessAssignment=
    FROM_SECTION)




if HolesSquare and ModelDim==3:
	r_4holes = 1.25
	d_4holes=7.5*0.8
	r_bearing = 5.5*0.75
	d_bearing = 4.

	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=1.8, name='__profile__', 
		sheetSize=72.0, transform=
		mdb.models['Model-1'].parts['UnitCellMain'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['UnitCellMain'].faces[16], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['UnitCellMain'].edges[6], 
		sketchOrientation=RIGHT, origin=(0.0, 0.0, 10.0)))
	mdb.models['Model-1'].parts['UnitCellMain'].projectReferencesOntoSketch(filter=
		COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
	mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(
		-11.7279220613579, 0.0), point2=(11.7279220613579, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
		addUndoState=False, entity=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[20], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[8], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
		11.7279220613579), point2=(0.0, -11.7279220613579))
	mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
		False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[21])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[25], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[21])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[14], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[21])
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		0.0, 4.95), point1=(2.25, 4.95))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[26], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[21])
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		6.3, 0.0), point1=(8.1, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[29], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[28], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		0.0, -5.85), point1=(1.8, -4.95))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[30], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[21])
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		-6.3, 0.0), point1=(-4.95, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[33], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[32], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].Spot(point=(0.0, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[34], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[20])
	mdb.models['Model-1'].sketches['__profile__'].EqualRadiusConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[25], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[22])
	mdb.models['Model-1'].sketches['__profile__'].EqualRadiusConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[22], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[23])
	mdb.models['Model-1'].sketches['__profile__'].EqualRadiusConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[23], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[24])
	mdb.models['Model-1'].sketches['__profile__'].EqualDistanceConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[26], entity2=
		mdb.models['Model-1'].sketches['__profile__'].vertices[28], midpoint=
		mdb.models['Model-1'].sketches['__profile__'].vertices[34])
	mdb.models['Model-1'].sketches['__profile__'].EqualDistanceConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[28], entity2=
		mdb.models['Model-1'].sketches['__profile__'].vertices[30], midpoint=
		mdb.models['Model-1'].sketches['__profile__'].vertices[34])
	mdb.models['Model-1'].sketches['__profile__'].EqualDistanceConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[30], entity2=
		mdb.models['Model-1'].sketches['__profile__'].vertices[32], midpoint=
		mdb.models['Model-1'].sketches['__profile__'].vertices[34])
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[34], entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[21])
	mdb.models['Model-1'].sketches['__profile__'].HorizontalDimension(textPoint=(
		3.9468412399292, -8.09041500091553), value=d_4holes, vertex1=
		mdb.models['Model-1'].sketches['__profile__'].vertices[34], vertex2=
		mdb.models['Model-1'].sketches['__profile__'].vertices[28])
	mdb.models['Model-1'].sketches['__profile__'].RadialDimension(curve=
		mdb.models['Model-1'].sketches['__profile__'].geometry[23], radius=r_4holes, 
		textPoint=(7.45872402191162, 6.54324436187744))
	mdb.models['Model-1'].parts['UnitCellMain'].CutExtrude(flipExtrudeDirection=OFF
		, sketch=mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=
		RIGHT, sketchPlane=mdb.models['Model-1'].parts['UnitCellMain'].faces[16], 
		sketchPlaneSide=SIDE1, sketchUpEdge=
		mdb.models['Model-1'].parts['UnitCellMain'].edges[6])
	del mdb.models['Model-1'].sketches['__profile__']
	
	
	

	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=1.83, name='__profile__', 
		sheetSize=73.23, transform=
		mdb.models['Model-1'].parts['UnitCellMain'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['UnitCellMain'].faces[5], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['UnitCellMain'].edges[37], 
		sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
	mdb.models['Model-1'].parts['UnitCellMain'].projectReferencesOntoSketch(filter=
		COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		0.0, 0.0), point1=(5.0325, 0.4575))
	mdb.models['Model-1'].sketches['__profile__'].RadialDimension(curve=
		mdb.models['Model-1'].sketches['__profile__'].geometry[23], radius=r_bearing, 
		textPoint=(6.13549423217773, 6.20962619781494))
	mdb.models['Model-1'].parts['UnitCellMain'].CutExtrude(depth=d_bearing, 
		flipExtrudeDirection=OFF, sketch=
		mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
		sketchPlane=mdb.models['Model-1'].parts['UnitCellMain'].faces[5], 
		sketchPlaneSide=SIDE1, sketchUpEdge=
		mdb.models['Model-1'].parts['UnitCellMain'].edges[37])
	del mdb.models['Model-1'].sketches['__profile__']
	mdb.models['Model-1'].parts['UnitCellMain'].Round(edgeList=(
		mdb.models['Model-1'].parts['UnitCellMain'].edges[1], ), radius=1.0)
#########################
#Part, Triangular cell
#########################
Partname = 'TriangularCell'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
#create corners

for ind in range(2):
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[1+3*ind], yv[1+3*ind]), point2=(
		xv[2+3*ind], yv[2+3*ind]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[1+3*ind], yv[1+3*ind]), point2=(
		xv[3+3*ind], yv[3+3*ind]))
	
	
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(xv[13], yv[13])
	, direction=COUNTERCLOCKWISE, point1=(xv[2], yv[2]), point2=(xv[6], yv[6]))


mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[17], yv[17]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[18], yv[18]))
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[17], yv[17]), point2=(
	xv[21], yv[21]))	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[18], yv[18]), point2=(
	xv[22], yv[22]))	
	
	
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[3], yv[3]), point2=(
	xv[21], yv[21]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[5], yv[5]), point2=(
	xv[22], yv[22]))


idP=1
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
    xv[idP]*0.999, yv[idP]*0.999), nearPoint2=(xv[idP]*0.999, 
    yv[idP]*0.999), radius=r_fillet)
	
idP=4
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[4], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[5], nearPoint1=(
    xv[idP]*0.999, yv[idP]*0.999), nearPoint2=(xv[idP]*0.999, 
    yv[idP]*0.999), radius=r_fillet)	
	
	
idP=0
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[7], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[8], nearPoint1=(
    0.0001, 0.0001), nearPoint2=(+0.0001, 
    +0.0001), radius=r_fillet)	
	
idP=21
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[9], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[11], nearPoint1=(
    xv[idP]*0.999, yv[idP]), nearPoint2=(xv[idP]*0.999, 
    yv[idP]), radius=r_fillet)	
		
idP=22
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[10], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[12], nearPoint1=(
    xv[idP], yv[idP]*0.999), nearPoint2=(xv[idP], 
    yv[idP]*0.999), radius=r_fillet)	
	



if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']



mdb.models['Model-1'].parts['TriangularCell'].PartitionEdgeByParam(edges=
    mdb.models['Model-1'].parts['TriangularCell'].edges[0:1]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[6:7]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[12:13], parameter=0.5)
	
	

mdb.models['Model-1'].parts['TriangularCell'].Set(edges=
    mdb.models['Model-1'].parts['TriangularCell'].edges[17:18], name=
    'SideEdgeTriangular')
mdb.models['Model-1'].parts['TriangularCell'].Set(edges=
    mdb.models['Model-1'].parts['TriangularCell'].edges[2:4]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[6:7]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[9:10]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[12:14]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[16:17]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[18:19], name=
    'TriangularLongStraightEdges')
mdb.models['Model-1'].parts['TriangularCell'].Set(edges=
    mdb.models['Model-1'].parts['TriangularCell'].edges[0:2]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[7:9]+\
    mdb.models['Model-1'].parts['TriangularCell'].edges[14:16], name=
    'TriangularHingeFilletEdges')


#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[0+dsec], thicknessAssignment=
    FROM_SECTION)




if MechanicalConnectionsTri and ModelDim ==3:

	R_circ = 0.9
	Edgedist1 = 1.6+1.7
	Edgedist2 = 1.4
	DoubleSpacing = -1.7
	flipdir=True
	depthcut=5.
	n1holes=7
	n2holes=n1holes-1
	Dist_circs = (Extrudedepth-2.*Edgedist2)/(n1holes-1)


	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.54, name='__profile__', 
		sheetSize=21.74, transform=
		mdb.models['Model-1'].parts['TriangularCell'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['TriangularCell'].faces[11], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['TriangularCell'].edges[35], 
		sketchOrientation=RIGHT, origin=(0., 0., 0.0)))
	mdb.models['Model-1'].parts['TriangularCell'].projectReferencesOntoSketch(
		filter=COPLANAR_EDGES, sketch=
		mdb.models['Model-1'].sketches['__profile__'])
		
	x1 = 	Edgedist1
	y1 =   Edgedist2
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		x1, y1), point1=(x1+R_circ, y1))
		
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		x1+DoubleSpacing, y1+Dist_circs/2), point1=(x1+DoubleSpacing+R_circ, y1+Dist_circs/2))		
		

			

	mdb.models['Model-1'].sketches['__profile__'].linearPattern(angle1=0.0, angle2=
		90.0, geomList=(mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((x1+R_circ, y1)), 
		), number1=1, number2=n1holes, spacing1=2.174, spacing2=Dist_circs, vertexList=())
		
	mdb.models['Model-1'].sketches['__profile__'].linearPattern(angle1=0.0, angle2=
		90.0, geomList=(mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((x1+DoubleSpacing+R_circ, y1+Dist_circs/2)), 
		), number1=1, number2=n2holes, spacing1=2.174, spacing2=Dist_circs, vertexList=())
			
	if 	flipdir:	
		mdb.models['Model-1'].parts['TriangularCell'].CutExtrude(depth=depthcut,flipExtrudeDirection=	
			ON, sketch=mdb.models['Model-1'].sketches['__profile__'], 
			sketchOrientation=RIGHT, sketchPlane=
			mdb.models['Model-1'].parts['TriangularCell'].faces[11], sketchPlaneSide=
			SIDE1, sketchUpEdge=
			mdb.models['Model-1'].parts['TriangularCell'].edges[35])
	
	else:	
		mdb.models['Model-1'].parts['TriangularCell'].CutExtrude(depth=depthcut,flipExtrudeDirection=	
			OFF, sketch=mdb.models['Model-1'].sketches['__profile__'], 
			sketchOrientation=RIGHT, sketchPlane=
			mdb.models['Model-1'].parts['TriangularCell'].faces[11], sketchPlaneSide=
			SIDE1, sketchUpEdge=
			mdb.models['Model-1'].parts['TriangularCell'].edges[35])

	del mdb.models['Model-1'].sketches['__profile__']




	





#########################
#Part, End Unit cell half cell
#########################


if ModelDim ==3:

	mdb.models['Model-1'].Part(name='UnitCellMain-RotatingHalf', objectToCopy=
		mdb.models['Model-1'].parts['UnitCellMain'])
	#del mdb.models['Model-1'].sketches['__profile__']
	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=2.02, name='__profile__', 
		sheetSize=81.01, transform=
		mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].faces.findAt(
		(0., dh3*0.9, Extrudedepth), ), sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].edges.findAt(
		(dh3, 1e-5, Extrudedepth), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 
		Extrudedepth)))
	"""	
	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=1.83, name='__profile__', 
		sheetSize=73.23, transform=
		mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].faces.findAt(
		(0., dh3, Extrudedepth), ), sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].edges.findAt(
		(0., dh3, Extrudedepth), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 
		Extrudedepth)))
	"""
	mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].projectReferencesOntoSketch(
		filter=COPLANAR_EDGES, sketch=
		mdb.models['Model-1'].sketches['__profile__'])
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
		-29.28, -11.4375))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-29.28, -11.4375), 
		point2=(1.3725, -27.45))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(1.3725, -27.45), 
		point2=(37.9725, 10.065))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(37.9725, 10.065), 
		point2=(8.6925, 26.9925))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(8.6925, 26.9925), 
		point2=(0.0, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
		0.0), point2=(-dh, dh))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((0.0, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((0.0, 0.0), )
		, entity2=mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((
		-6.363961, 6.363961), ))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((-dh, 
		dh))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
		addUndoState=False, entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((
		-dh, dh), ), entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961), ))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-14.64, 
		-5.71875))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].AngularDimension(line1=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-14.64, 
		-5.71875), ), line2=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961), ), textPoint=(-13.3068599700928, 5.04662322998047), value=60.0)
	mdb.models['Model-1'].sketches['__profile__'].undo()
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-14.64, 
		-5.71875))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].AngularDimension(line1=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-14.64, 
		-5.71875), ), line2=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961), ), textPoint=(-15.2351055145264, 6.77726078033447), value=65.0)
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((-30.160515, 
		-10.97753))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((8.6925, 
		26.9925))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].SymmetryConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((
		-30.1605147376598, -10.9775296146553), ), entity2=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((8.6925, 
		26.9925), ), symmetryAxis=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961), ))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((1.3725, -27.45))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((37.9725, 
		10.065))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].SymmetryConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((1.3725, 
		-27.45), ), entity2=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((37.9725, 
		10.065), ), symmetryAxis=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961), ))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((0.0, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((-30.160515, 
		-10.97753))
	mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
		-32.3251647949219, 7.20425510406494), value=50.0, vertex1=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((0.0, 0.0), )
		, vertex2=mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((
		-30.1605147376598, -10.9775296146553), ))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-28.524816, 
		-27.536754))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961))
	mdb.models['Model-1'].sketches['__profile__'].ParallelConstraint(entity1=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-28.524816, 
		-27.536754), ), entity2=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-6.363961, 
		6.363961), ))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((0.0, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((0.0, 0.0), 
		))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((-46.984631, 
		-17.101007))
	mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((-16.995647, 
		-47.089992))
	mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
		-54.3123588562012, -29.9055309295654), value=50.0, vertex1=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((
		-46.9846310392954, -17.1010071662835), ), vertex2=
		mdb.models['Model-1'].sketches['__profile__'].vertices.findAt((
		-16.9956466518108, -47.089991553768), ))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-23.492316, 
		-8.550504))
	mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((8.550504, 
		23.492316))
	mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((-23.492316, 
		-8.550504), ), curve2=
		mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((8.550504, 
		23.492316), ), nearPoint1=(-4.24820899963379, -2.10990285873413), 
		nearPoint2=(1.18997383117676, 4.21935653686523), radius=1.0)
	mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].CutExtrude(
		flipExtrudeDirection=OFF, sketch=
		mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
		sketchPlane=
		mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].faces.findAt((
		1e-5, 0.9*dh3, Extrudedepth), ), sketchPlaneSide=SIDE1, sketchUpEdge=
		mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].edges.findAt((
		dh3, 1e-5, Extrudedepth), ))
	del mdb.models['Model-1'].sketches['__profile__']




	#Section assignment
	mdb.models['Model-1'].parts['UnitCellMain-RotatingHalf'].SectionAssignment(offset=0.0, 
		offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
		faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
		mask=('[#1 ]', ), )), sectionName=SectionsAll[0+dsec], thicknessAssignment=
		FROM_SECTION)










































# Create Triangular Hinge Geometry
Partname = 'HingeTriangularCell'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[25+10], yv[25+10]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[26+10], yv[26+10]))
	

mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
    0.0001, 0.0001), nearPoint2=(+0.0001, 
    +0.0001), radius=r_fillet_tri)	
	


mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[27+10], yv[27+10]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[28+10], yv[28+10]))
	
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[5], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[6], nearPoint1=(
    -0.0001, -0.0001), nearPoint2=(-0.0001, 
    -0.0001), radius=r_fillet_tri)	
		
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[25+10], yv[25+10]), point2=(
	xv[28+10], yv[28+10]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26+10], yv[26+10]), point2=(
	xv[27+10], yv[27+10]))
	
	

if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']

if CutwTriangularHinge:

	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.74, name='__profile__', 
		sheetSize=29.74, transform=
		mdb.models['Model-1'].parts['HingeTriangularCell'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['HingeTriangularCell'].faces[8], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['HingeTriangularCell'].edges[16], 
		sketchOrientation=RIGHT, origin=(0.0, 0.0, 10.0)))
	mdb.models['Model-1'].parts['HingeTriangularCell'].projectReferencesOntoSketch(
		filter=COPLANAR_EDGES, sketch=
		mdb.models['Model-1'].sketches['__profile__'])
	mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-1000.,1000.), 
		point2=(-w_hinge_tri/2, -1000.))
	mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(1000.,1000.), 
		point2=(w_hinge_tri/2, -1000.))
	mdb.models['Model-1'].parts['HingeTriangularCell'].CutExtrude(
		flipExtrudeDirection=OFF, sketch=
		mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
		sketchPlane=mdb.models['Model-1'].parts['HingeTriangularCell'].faces[8], 
		sketchPlaneSide=SIDE1, sketchUpEdge=
		mdb.models['Model-1'].parts['HingeTriangularCell'].edges[16])
	del mdb.models['Model-1'].sketches['__profile__']

#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[2+dsec], thicknessAssignment=
    FROM_SECTION)
	
	
	
	
	





# Create Square Hinge Geometry
Partname = 'HingeSquareVertical'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
if not CutHingesSquare:
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[25], yv[25]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[26], yv[26]))
		

	mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
		0.0001, 0.0001), nearPoint2=(+0.0001, 
		+0.0001), radius=r_fillet)	
		


	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[27], yv[27]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[28], yv[28]))
		
	mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[5], curve2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[6], nearPoint1=(
		-0.0001, -0.0001), nearPoint2=(-0.0001, 
		-0.0001), radius=r_fillet)	
			
		
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[25], yv[25]), point2=(
		xv[28], yv[28]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26], yv[26]), point2=(
		xv[27], yv[27]))

	mdb.models['Model-1'].sketches['__profile__'].rotate(angle=45.0, centerPoint=(0.0, 
		0.0), objectList=maketuple(mdb.models['Model-1'].sketches['__profile__'].geometry))

elif CutHingesSquare:
	dy = pv2*Stretch_fac_HingesSq
	dx = w_hinge*wrat_hinge
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-dx/2, dy/2), point2=(dx/2,dy/2))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(dx/2, dy/2), point2=(dx/2,-dy/2))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(dx/2, -dy/2), point2=(-dx/2,-dy/2))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-dx/2, -dy/2), point2=(-dx/2,dy/2))



if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']


#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[1+dsec], thicknessAssignment=
    FROM_SECTION)	








# Create Square Hinge Geometry
Partname = 'HingeSquareHorizontal'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
if not CutHingesSquare:

	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[25], yv[25]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[26], yv[26]))
		

	mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
		0.0001, 0.0001), nearPoint2=(+0.0001, 
		+0.0001), radius=r_fillet)	
		


	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[27], yv[27]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
		xv[28], yv[28]))
		
	mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
		mdb.models['Model-1'].sketches['__profile__'].geometry[5], curve2=
		mdb.models['Model-1'].sketches['__profile__'].geometry[6], nearPoint1=(
		-0.0001, -0.0001), nearPoint2=(-0.0001, 
		-0.0001), radius=r_fillet)	
			
		
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[25], yv[25]), point2=(
		xv[28], yv[28]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26], yv[26]), point2=(
		xv[27], yv[27]))
	mdb.models['Model-1'].sketches['__profile__'].rotate(angle=-45.0, centerPoint=(0.0, 
		0.0), objectList=maketuple(mdb.models['Model-1'].sketches['__profile__'].geometry))

elif CutHingesSquare:
	dx = pv2*Stretch_fac_HingesSq
	dy = w_hinge*wrat_hinge
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-dx/2, dy/2), point2=(dx/2,dy/2))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(dx/2, dy/2), point2=(dx/2,-dy/2))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(dx/2, -dy/2), point2=(-dx/2,-dy/2))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-dx/2, -dy/2), point2=(-dx/2,dy/2))




if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']


#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[1+dsec], thicknessAssignment=
    FROM_SECTION)	



# Create Square Hinge Mold part 1Geometry
Partname = 'MOLD1HingeSquareHorizontal'

t_mold=2.


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[25], yv[25]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[26], yv[26]))
	

mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
    0.0001, 0.0001), nearPoint2=(+0.0001, 
    +0.0001), radius=r_fillet)	
	

"""
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[27], yv[27]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[28], yv[28]))
	
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[5], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[6], nearPoint1=(
    -0.0001, -0.0001), nearPoint2=(-0.0001, 
    -0.0001), radius=r_fillet)	
"""
		
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[25], yv[25]), point2=(
	(xv[28]+xv[25])/2, (yv[28]+yv[25])/2))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26], yv[26]), point2=(
	(xv[27]+xv[26])/2, (yv[27]+yv[26])/2))
"""	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(
	(xv[28]+xv[25])/2, (yv[28]+yv[25])/2), point2=(
	(xv[27]+xv[26])/2, (yv[27]+xv[27])/2))	
"""




mdb.models['Model-1'].sketches['__profile__'].rotate(angle=-45.0, centerPoint=(0.0,   0.0), objectList=maketuple(mdb.models['Model-1'].sketches['__profile__'].geometry))
	
	
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,wh2), point2=(-t_mold/2.,wh2+t_mold/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,wh2+t_mold), point2=(-t_mold/2.,wh2+t_mold/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,wh2+t_mold), point2=(lh2+t_mold,wh2+t_mold))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,wh2), point2=(lh2+t_mold,wh2+t_mold))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,wh2), point2=(lh2+t_mold-wh2,0.))

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,-wh2), point2=(lh2+t_mold-wh2,0.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,-wh2), point2=(lh2+t_mold,-(wh2+t_mold)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,-(wh2+t_mold)), point2=(lh2+t_mold,-(wh2+t_mold)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,-(wh2+t_mold)), point2=(-t_mold/2.,-(wh2+t_mold/2.)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,-wh2), point2=(-t_mold/2.,-(wh2+t_mold/2.)))
	



if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']


#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[1+dsec], thicknessAssignment=
    FROM_SECTION)	





# Create Square Hinge Mold part 2 Geometry
Partname = 'MOLD2HingeSquareHorizontal'

t_mold=2.


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[25], yv[25]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[26], yv[26]))
	

mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], nearPoint1=(
    0.0001, 0.0001), nearPoint2=(+0.0001, 
    +0.0001), radius=r_fillet)	
	

"""
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[27], yv[27]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[0], yv[0]), point2=(
	xv[28], yv[28]))
	
mdb.models['Model-1'].sketches['__profile__'].FilletByRadius(curve1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[5], curve2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[6], nearPoint1=(
    -0.0001, -0.0001), nearPoint2=(-0.0001, 
    -0.0001), radius=r_fillet)	
"""
		
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[25], yv[25]), point2=(
	(xv[28]+xv[25])/2, (yv[28]+yv[25])/2))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26], yv[26]), point2=(
	(xv[27]+xv[26])/2, (yv[27]+yv[26])/2))
"""	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(
	(xv[28]+xv[25])/2, (yv[28]+yv[25])/2), point2=(
	(xv[27]+xv[26])/2, (yv[27]+xv[27])/2))	
"""




mdb.models['Model-1'].sketches['__profile__'].rotate(angle=-45.0, centerPoint=(0.0,   0.0), objectList=maketuple(mdb.models['Model-1'].sketches['__profile__'].geometry))
	
	
	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,wh2), point2=(t_mold/2.,wh2+t_mold/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,wh2+t_mold), point2=(t_mold/2.,wh2+t_mold/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,wh2+t_mold), point2=(lh2+t_mold,wh2+t_mold))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,wh2), point2=(lh2+t_mold,wh2+t_mold))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,wh2), point2=(lh2+t_mold+wh2,0.))

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,-wh2), point2=(lh2+t_mold+wh2,0.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(lh2+t_mold,-wh2), point2=(lh2+t_mold,-(wh2+t_mold)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,-(wh2+t_mold)), point2=(lh2+t_mold,-(wh2+t_mold)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,-(wh2+t_mold)), point2=(t_mold/2.,-(wh2+t_mold/2.)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.,-wh2), point2=(t_mold/2.,-(wh2+t_mold/2.)))
	



if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']


#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[1+dsec], thicknessAssignment=
    FROM_SECTION)	




##########################################
# Create ConnectionHinge for Casting
##########################################
if TypeHingeTriangle==1:
	Partname = 'ConnectionCasting'


	mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26], yv[26]), point2=(
		xv[27], yv[27]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[26], yv[26]), point2=(
		xv[30], yv[30]))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[27], yv[27]), point2=(
		xv[31], yv[31]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[30], yv[30]), point2=(
		xv[33], yv[33]))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(xv[31], yv[31]), point2=(
		xv[33], yv[33]))
		



	if ModelDim==2:
		mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
			, type=DEFORMABLE_BODY)
		mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
			mdb.models['Model-1'].sketches['__profile__'])
	elif ModelDim == 3:
		mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
			, type=DEFORMABLE_BODY)	
		mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=Extrudedepth, sketch=
			mdb.models['Model-1'].sketches['__profile__'])		
	del mdb.models['Model-1'].sketches['__profile__']


	#Section assignment
	mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
		offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
		faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
		mask=('[#1 ]', ), )), sectionName=SectionsAll[0+dsec], thicknessAssignment=
		FROM_SECTION)
		
		


##########################################
# Create MoldSupport
##########################################

Partname = 'MoldSupport'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

mdb.models['Model-1'].sketches['__profile__'].Line(point1=((wv1/2.+2.*rv1), (wv1/2.)), point2=(
	(wv1/2.+2.*rv1), -(wv1/2.)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=((wv1/2.+2.*rv1), -(wv1/2.)), point2=(
	-(wv1/2.+2.*rv1), -(wv1/2.)))	
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-(wv1/2.+2.*rv1), -(wv1/2.)), point2=(
	-(wv1/2.+2.*rv1), (wv1/2.)))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-(wv1/2.+2.*rv1), (wv1/2.)), point2=(
	(wv1/2.+2.*rv1), (wv1/2.)))
	



if ModelDim==2:
	mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=Partname
		, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[Partname].BaseShell(sketch=
		mdb.models['Model-1'].sketches['__profile__'])
elif ModelDim == 3:
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Partname
		, type=DEFORMABLE_BODY)	
	mdb.models['Model-1'].parts[Partname].BaseSolidExtrude(depth=depthv1, sketch=
		mdb.models['Model-1'].sketches['__profile__'])		
del mdb.models['Model-1'].sketches['__profile__']
if ModelDim == 3:
	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=2.39, name='__profile__', 
		sheetSize=95.69, transform=
		mdb.models['Model-1'].parts[Partname].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts[Partname].faces[4], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts[Partname].edges[7], 
		sketchOrientation=RIGHT, origin=(0.0, 0.0,depthv1 )))
	mdb.models['Model-1'].parts[Partname].projectReferencesOntoSketch(filter=
		COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
		
		
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-pv2/2, wh2), point2=((pv2/2, wh2)))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-pv2/2, wh2), point2=((-pv2/2, wh2+0.1*tv1)))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(pv2/2, wh2), point2=((pv2/2, wh2+0.1*tv1)))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-pv2/2, wh2+0.1*tv1), point2=((0., wh2+tv1))	)
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(pv2/2, wh2+0.1*tv1), point2=((0., wh2+tv1))	)
	
		
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-pv2/2, -wh2), point2=((pv2/2, -wh2)))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-pv2/2, -wh2), point2=((-pv2/2, -(wh2+0.1*tv1))))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(pv2/2, -wh2), point2=((pv2/2, -(wh2+0.1*tv1))))	
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-pv2/2, -(wh2+0.1*tv1)), point2=((0., -(wh2+tv1)))	)
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(pv2/2, -(wh2+0.1*tv1)), point2=((0., -(wh2+tv1)))	)
	"""	
	mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-pv2/2, wh2), 
		point2=(pv2/2, wh2+tv1))	
	mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-pv2/2, -wh2), 
		point2=(pv2/2, -(wh2+tv1)))
	"""
		
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		-(wv1/2.+rv1), 0.0), point1=(-wv1/2., 0.))
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		(wv1/2.+rv1), 0.0), point1=(wv1/2., 0.))

	mdb.models['Model-1'].parts[Partname].SolidExtrude(depth=Extrudedepth, 
		flipExtrudeDirection=OFF, sketch=
		mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
		sketchPlane=mdb.models['Model-1'].parts[Partname].faces[4], 
		sketchPlaneSide=SIDE1, sketchUpEdge=
		mdb.models['Model-1'].parts[Partname].edges[7])
	del mdb.models['Model-1'].sketches['__profile__']

#Section assignment
mdb.models['Model-1'].parts[Partname].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    faces=mdb.models['Model-1'].parts[Partname].faces.getSequenceFromMask(
    mask=('[#1 ]', ), )), sectionName=SectionsAll[0+dsec], thicknessAssignment=
    FROM_SECTION)
	
	


##########################################
# Create Mold Assembly
##########################################
if CreateMoldAssembly:
	mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
		'TriangularCell-1', part=mdb.models['Model-1'].parts['TriangularCell'])
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
		'TriangularCell-2', part=mdb.models['Model-1'].parts['TriangularCell'])
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='MoldSupport-1', 
		part=mdb.models['Model-1'].parts['MoldSupport'])
	mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, 
		1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('TriangularCell-2', ))
	mdb.models['Model-1'].rootAssembly.rotate(angle=45.0, axisDirection=(0.0, 0.0, 
		1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('MoldSupport-1', ))
	mdb.models['Model-1'].rootAssembly.translate(instanceList=('MoldSupport-1', ), 
		vector=(0.0, 0.0, -depthv1))
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=(
		mdb.models['Model-1'].rootAssembly.instances['TriangularCell-1'], 
		mdb.models['Model-1'].rootAssembly.instances['TriangularCell-2'], 
		mdb.models['Model-1'].rootAssembly.instances['MoldSupport-1']), 
		keepIntersections=ON, name='MoldAssembly', originalInstances=DELETE)
	del mdb.models['Model-1'].rootAssembly.features['MoldAssembly-1']


##########################################
# Create Rigid Top
##########################################


PartName = 'Rigid_Top'


mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-d_cell, (Ny-1)*d_cell+dh3+h_rigid), point2=    (Nx*d_cell, (Ny-1)*d_cell+dh3+h_rigid))

mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name=PartName, type=
    DISCRETE_RIGID_SURFACE)
mdb.models['Model-1'].parts[PartName].BaseWire(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].parts[PartName].ReferencePoint(point=((Nx-1)/2*d_cell, (Ny-1)*d_cell+dh3+h_rigid+h_rigid_RP, 0.0))
mdb.models['Model-1'].parts['Rigid_Top'].Set(name='RP_Rigid', referencePoints=(
    mdb.models['Model-1'].parts[PartName].referencePoints[mdb.models['Model-1'].parts[PartName].referencePoints.keys()[0]], ))
if Solvetype=='Implicit':	
	mdb.models['Model-1'].parts['Rigid_Top'].setElementType(elemTypes=(ElemType(
		elemCode=R2D2, elemLibrary=STANDARD), ), regions=(
		mdb.models['Model-1'].parts['Rigid_Top'].edges.findAt(((0., (Ny-1)*d_cell+dh3+h_rigid, 0.0), )), ))
elif Solvetype=='Explicit':	
	mdb.models['Model-1'].parts['Rigid_Top'].setElementType(elemTypes=(ElemType(
		elemCode=R2D2, elemLibrary=EXPLICIT), ), regions=(
		mdb.models['Model-1'].parts['Rigid_Top'].edges.findAt(((0., (Ny-1)*d_cell+dh3+h_rigid, 0.0), )), ))


	
mdb.models['Model-1'].parts[PartName].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=MeshSize)
mdb.models['Model-1'].parts[PartName].generateMesh()


mdb.models['Model-1'].parts['Rigid_Top'].engineeringFeatures.PointMassInertia(
    alpha=alpha_contact, composite=0.0, mass=RigidMass/2, name='Mass', region=Region(
    vertices=mdb.models['Model-1'].parts['Rigid_Top'].vertices.findAt(((-d_cell, (Ny-1)*d_cell+dh3+h_rigid,0.),), ((Nx*d_cell, (Ny-1)*d_cell+dh3+h_rigid,0.), ), )))


#####################################
# Create Merged Parts
####################################
PartName1 = 'TriangularCell'
PartName2 = 'TriangularCell'
if TypeHingeTriangle==0:	
	PartName3 = 'HingeTriangularCell'
elif TypeHingeTriangle==1:	
	PartName3 = 'ConnectionCasting'
	PartName4 = 'ConnectionCasting'

InsName1 = PartName1+'-1'
InsName2 = PartName2+'-2'
if TypeHingeTriangle==0:	
	InsName3 = PartName3+'-1'
elif TypeHingeTriangle==1:	
	InsName3 = PartName3+'-1'
	InsName4 = PartName4+'-2'


CombPartName = 'DoubleTriangularCellHinged'
CombInsName = CombPartName+'-1'

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
    InsName1, part=
    mdb.models['Model-1'].parts[PartName1])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
    InsName2, part=
    mdb.models['Model-1'].parts[PartName2])
if TypeHingeTriangle==0 or TypeHingeTriangle==1:	
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
		InsName3, part=
		mdb.models['Model-1'].parts[PartName3])	


mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, 
    1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=(InsName2, ))



if TypeHingeTriangle==1:	
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
		InsName4, part=
		mdb.models['Model-1'].parts[PartName4])		
	mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, 
		1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=(InsName4, ))

if TypeHingeTriangle==0:		
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=(
		mdb.models['Model-1'].rootAssembly.instances[InsName1],
		mdb.models['Model-1'].rootAssembly.instances[InsName2], 
		mdb.models['Model-1'].rootAssembly.instances[InsName3]), 
		keepIntersections=ON, name=CombPartName, originalInstances=
		DELETE)
elif TypeHingeTriangle==1:	
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=(
		mdb.models['Model-1'].rootAssembly.instances[InsName1],
		mdb.models['Model-1'].rootAssembly.instances[InsName2], 
		mdb.models['Model-1'].rootAssembly.instances[InsName3], 
		mdb.models['Model-1'].rootAssembly.instances[InsName4]), 
		keepIntersections=OFF, name=CombPartName, originalInstances=
		DELETE)
elif TypeHingeTriangle==2:	
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=(
		mdb.models['Model-1'].rootAssembly.instances[InsName1],
		mdb.models['Model-1'].rootAssembly.instances[InsName2]), 
		keepIntersections=OFF, name=CombPartName, originalInstances=
		DELETE)
	
del mdb.models['Model-1'].rootAssembly.features[CombInsName]



#####################################
# Create Assembly
####################################	

PartListMain = ['UnitCellMain', CombPartName,'HingeTriangularCell','UnitCellMain-RotatingHalf']
InsNamestart = 'Cell'

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

for indX in range(Nx):
	for indY in range(Ny):
		if Assembly_mat[indX,indY] > -1e-10:
			PartName = PartListMain[Assembly_mat[indX,indY]]
			InsName = InsNamestart+'-'+str(indX)+'-'+str(indY)
			TX = indX*d_cell
			TY = indY*d_cell
			mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
				InsName, part=    mdb.models['Model-1'].parts[PartName])
			mdb.models['Model-1'].rootAssembly.rotate(angle=Rotate_mat[indX,indY], axisDirection=(0.0, 0.0, 1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=(InsName, ))			
			mdb.models['Model-1'].rootAssembly.rotate(angle=Rotate_mat_mainonly[indX,indY], axisDirection=(0.0, 0.0, 1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=(InsName, ))	
			
			mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))


if MergeIndividualParts:
	AssemblyName = 'MergeMainParts'
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=maketuple(mdb.models['Model-1'].rootAssembly.instances), 
		keepIntersections=ON, name=AssemblyName, originalInstances=
		DELETE)
	del mdb.models['Model-1'].rootAssembly.features[AssemblyName+'-1']



if AssembleTriangularHinges:
	PartName = PartListMain[2]
	InsNamestart = 'TriangularHinge'
	for indX in range(Nx):
		for indY in range(Ny):
			if Assembly_mat[indX,indY]==1:
				InsName = InsNamestart+'-'+str(indX)+'-'+str(indY)
				TX = indX*d_cell
				TY = indY*d_cell
				mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
					InsName, part=    mdb.models['Model-1'].parts[PartName])							
				mdb.models['Model-1'].rootAssembly.rotate(angle=Rotate_mat[indX,indY], axisDirection=(0.0, 0.0, 1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=(InsName, ))
				mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))




	AssemblyName = 'AssemblyTriangularHinges'
	if len(mdb.models['Model-1'].rootAssembly.instances) >1.1:	
		mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
			instances=maketuple(mdb.models['Model-1'].rootAssembly.instances), 
			keepIntersections=ON, name=AssemblyName, originalInstances=
			DELETE)
		del mdb.models['Model-1'].rootAssembly.features[AssemblyName+'-1']
				



			
if TypeHingeTriangle==1 and MergeAssembly == False:
	hfsay		
if AssembleHinges:
	PartName = 'HingeSquareVertical'
	InsNamestart = 'HingeVertical'
	for indX in range(Nx):
		for indY in range(Ny-1):
			if Assembly_mat[indX,indY] > -1e-10 and Assembly_mat[indX,indY+1] > -1e-10: 
				InsName = InsNamestart+'-'+str(indX)+'-'+str(indY)
				TX = indX*d_cell
				TY = (indY+0.5)*(d_cell)
				mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
					InsName, part=    mdb.models['Model-1'].parts[PartName])
							
				mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))
	if CutHingesSquare:
		Disp_mat_Vertical =-np.tan(np.deg2rad(Rotate_mat))	* d_cell/2			
		for indX in range(Nx):
			for indY in range(Ny-1):
				if Assembly_mat[indX,indY] > -1e-10 and Assembly_mat[indX,indY+1] > -1e-10: 	
						InsName = InsNamestart+'-'+str(indX)+'-'+str(indY)
						TX = Disp_mat_Vertical[indX,indY]
						TY = 0.
						mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))				
					
				



	
	PartName = 'HingeSquareHorizontal'
	InsNamestart = 'HingeHorizontal'
		
	for indX in range(Nx-1):
		for indY in range(Ny):
			if Assembly_mat[indX,indY] > -1e-10 and Assembly_mat[indX+1,indY] > -1e-10: 
				InsName = InsNamestart+'-'+str(indX)+'-'+str(indY)
				TX = (indX+0.5)*(d_cell)
				TY = (indY)*(d_cell)
				mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
					InsName, part=    mdb.models['Model-1'].parts[PartName])
							
				mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))		
	if CutHingesSquare:
		Disp_mat_Horizontal =np.tan(np.deg2rad(Rotate_mat))	* d_cell/2			
		for indX in range(Nx-1):
			for indY in range(Ny):
				if Assembly_mat[indX,indY] > -1e-10 and Assembly_mat[indX+1,indY] > -1e-10: 	
						InsName = InsNamestart+'-'+str(indX)+'-'+str(indY)
						TX = 0.
						TY = Disp_mat_Horizontal[indX,indY]
						mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))		


if MergeIndividualParts:
	AssemblyName = 'MergeSquareHinges'
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=maketuple(mdb.models['Model-1'].rootAssembly.instances), 
		keepIntersections=ON, name=AssemblyName, originalInstances=
		DELETE)
	del mdb.models['Model-1'].rootAssembly.features[AssemblyName+'-1']






































if MergeIndividualParts:	
	PartName = 'MergeMainParts'
	InsName = PartName+'-1'
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
		InsName, part=    mdb.models['Model-1'].parts[PartName])	

	PartName = 'MergeSquareHinges'
	InsName = PartName+'-1'
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
		InsName, part=    mdb.models['Model-1'].parts[PartName])


	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
		mdb.models['Model-1'].rootAssembly.instances['MergeMainParts-1'], ), 
		instanceToBeCut=
		mdb.models['Model-1'].rootAssembly.instances['MergeSquareHinges-1'], name=
		'SquareHingesCut', originalInstances=DELETE)

	PartName = 'MergeMainParts'
	InsName = PartName+'-1'
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=	
		InsName, part=    mdb.models['Model-1'].parts[PartName])
















































AssemblyName = 'AssemblyPart'
if len(mdb.models['Model-1'].rootAssembly.instances) >1.1:	
	mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
		instances=maketuple(mdb.models['Model-1'].rootAssembly.instances), 
		keepIntersections=ON, name=AssemblyName, originalInstances=
		DELETE)
	#del mdb.models['Model-1'].rootAssembly.features[AssemblyName+'-1']
			



if RotatedOrientation:
	mdb.models['Model-1'].rootAssembly.regenerate()
	mdb.models['Model-1'].rootAssembly.rotate(angle=-45.0, axisDirection=(0.0, 0.0, 
		1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('AssemblyPart-1', ))

##################CONTINUE HERE


# ADD Rigid 
PartName = 'Rigid_Top'
if BCTop == 'rigid':
	InsName = 'Rigid-1'
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=InsName, 
		part=mdb.models['Model-1'].parts[PartName])

if BCBottom == 'rigid':
	InsName = 'Rigid-2'
	TX = 0
	TY = -(Ny-1)*d_cell-2*dh3-2.*h_rigid
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=InsName, 
		part=mdb.models['Model-1'].parts[PartName])	
	mdb.models['Model-1'].rootAssembly.translate(instanceList=(InsName, ), vector=(TX, TY, 0.0))		





if WriteSTP:
	try:
		part_name = 'MergeMainParts'
		stp_name = file_location+'\\Model_v0'+str(Model_version)+'_0'+str(Subversion)+'_'+part_name+'.stp'
		mdb.models['Model-1'].parts[part_name].writeStepFile(stp_name)

		part_name = 'AssemblyTriangularHinges'
		stp_name = file_location+'\\Model_v0'+str(Model_version)+'_0'+str(Subversion)+'_'+part_name+'.stp'
		mdb.models['Model-1'].parts[part_name].writeStepFile(stp_name)
		
		part_name = 'SquareHingesCut'
		stp_name = file_location+'\\Model_v0'+str(Model_version)+'_0'+str(Subversion)+'_'+part_name+'.stp'
		mdb.models['Model-1'].parts[part_name].writeStepFile(stp_name)		


		mdb.saveAs(pathName=file_name)	
	except:
		pass


###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
#####################################
# NOTE: THE FOLLOWING DEALS WITH THE ANALYSIS, NOT THE GEOMETRY
####################################	
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################













if TypeHingeTriangle==0:	
	
		
		
	#####################################
	# Mesh
	####################################		
	PartName = 'AssemblyPart'
	mdb.models['Model-1'].parts['AssemblyPart'].seedPart(  size=MeshSize)

	#mdb.models['Model-1'].parts['AssemblyPart'].seedPart(deviationFactor=0.02, 
	#	minSizeFactor=0.02, size=MeshSize)
	if SeedEdge:
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(edges=
			mdb.models['Model-1'].parts['AssemblyPart'].sets['SideEdges'].edges[:], number=SeedCountSide)
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(edges=
			mdb.models['Model-1'].parts['AssemblyPart'].sets['SideEdgeTriangular'].edges[:], number=SeedCountSide)
			
		"""
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(edges=
			mdb.models['Model-1'].parts['AssemblyPart'].sets['CornerStraightEdges'].edges[:], number=SeedCountStraight)
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(edges=
			mdb.models['Model-1'].parts['AssemblyPart'].sets['TriangularLongStraightEdges'].edges[:], number=SeedCountStraight)
		"""
			
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(edges=
			mdb.models['Model-1'].parts['AssemblyPart'].sets['HingeEdges'].edges[:], number=SeedCountHinge)
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(edges=
			mdb.models['Model-1'].parts['AssemblyPart'].sets['TriangularHingeFilletEdges'].edges[:], number=SeedCountHinge)
			
		#mdb.models['Model-1'].parts['AssemblyPart'].seedEdgeByNumber(edges=
		#	mdb.models['Model-1'].parts['AssemblyPart'].set['SideEdges'].edges[:], number=SeedCount)
			
		"""
		mdb.models['Model-1'].parts[PartName].seedEdgeByNumber(constraint=FINER, edges=mdb.models['Model-1'].parts[PartName].edges[:], number=SeedCount)
		"""
	

	
	ElType = 'Quadratic'
	if Solvetype =='Implicit':
		if ElType == 'Linear':	
			mdb.models['Model-1'].parts['AssemblyPart'].setElementType(elemTypes=(ElemType(
				elemCode=CPE4, elemLibrary=STANDARD), ElemType(elemCode=CPE3, 
				elemLibrary=STANDARD)), regions=(
				mdb.models['Model-1'].parts['AssemblyPart'].faces[:], ))
		elif ElType == 'Quadratic':	
			mdb.models['Model-1'].parts['AssemblyPart'].setMeshControls(elemShape=TRI, 
				regions=mdb.models['Model-1'].parts['AssemblyPart'].faces[:])
			mdb.models['Model-1'].parts['AssemblyPart'].setElementType(elemTypes=(ElemType(
				elemCode=CPE6, elemLibrary=STANDARD), ElemType(elemCode=CPE6, 
				elemLibrary=STANDARD)), regions=(mdb.models['Model-1'].parts['AssemblyPart'].faces[:], ))
	elif Solvetype =='Explicit':
		if ElType == 'Linear':	
			mdb.models['Model-1'].parts['AssemblyPart'].setElementType(elemTypes=(ElemType(
				elemCode=CPE4R, elemLibrary=EXPLICIT), ElemType(elemCode=CPE3, 
				elemLibrary=EXPLICIT)), regions=(
				mdb.models['Model-1'].parts['AssemblyPart'].faces[:], ))
		elif ElType == 'Quadratic':	
			mdb.models['Model-1'].parts['AssemblyPart'].setMeshControls(elemShape=
				TRI, regions=
				mdb.models['Model-1'].parts['AssemblyPart'].faces[:])
			mdb.models['Model-1'].parts['AssemblyPart'].setElementType(elemTypes=(ElemType(
				elemCode=UNKNOWN_QUAD, elemLibrary=EXPLICIT), ElemType(elemCode=CPE6M, 
				elemLibrary=EXPLICIT)), regions=(
				mdb.models['Model-1'].parts['AssemblyPart'].faces[:], ))	
	
	mdb.models['Model-1'].parts['AssemblyPart'].generateMesh()




		
	#####################################
	# Define Sets
	####################################

	Points =[]
	for indX in range(Nx):
		TX = (indX)*(d_cell)
		TY = -dh3
		Points = Points+[mdb.models['Model-1'].parts['AssemblyPart'].vertices.findAt(((TX,TY 
		, 0.0), ))]
	Points = tuple(Points)
	mdb.models['Model-1'].parts['AssemblyPart'].Set(name='Bottom', vertices=Points)




	Points =[]
	for indX in range(Nx):
		TX = (indX)*(d_cell)
		TY = (Ny-1)*d_cell+dh3
		Points = Points+[mdb.models['Model-1'].parts['AssemblyPart'].vertices.findAt(((TX,TY 
		, 0.0), ))]
	Points = tuple(Points)
	mdb.models['Model-1'].parts['AssemblyPart'].Set(name='Top', vertices=Points)




	Points =[]
	for indY in range(Ny):
		TX = -dh3
		TY = (indY)*(d_cell)
		Points = Points+[mdb.models['Model-1'].parts['AssemblyPart'].vertices.findAt(((TX,TY 
		, 0.0), ))]
	Points = tuple(Points)
	mdb.models['Model-1'].parts['AssemblyPart'].Set(name='Left', vertices=Points)



	Points =[]
	for indY in range(Ny):
		TX = (Nx-1)*d_cell+dh3
		TY = (indY)*(d_cell)
		Points = Points+[mdb.models['Model-1'].parts['AssemblyPart'].vertices.findAt(((TX,TY 
		, 0.0), ))]
	Points = tuple(Points)
	mdb.models['Model-1'].parts['AssemblyPart'].Set(name='Right', vertices=Points)


		
	#####################################
	# Define Step
	####################################
	###########################
	# Step
	###########################	
		
	if Solvetype=='Implicit':
		if steptype =='Buckling':
			mdb.models['Model-1'].BuckleStep(blockSize=DEFAULT, eigensolver=LANCZOS, 
				maxBlocks=DEFAULT, minEigen=0, name='Loading', numEigen=10, 
				previous='Initial')
			
		elif steptype =='StaticNL':
				
			mdb.models['Model-1'].StaticStep(adaptiveDampingRatio=None, 
				continueDampingFactors=False, initialInc=1./MinSteps, matrixSolver=DIRECT, 
				matrixStorage=SOLVER_DEFAULT, maxInc=1./MinSteps, maxNumInc=1000, minInc=1e-05/MinSteps, 
				name='Loading', nlgeom=ON, previous='Initial', solutionTechnique=
				FULL_NEWTON, stabilizationMethod=NONE)
		elif steptype =='ViscoStep':
			t_loading = t_loading_vec[0]
			maxincloading = t_loading/MinSteps 
			mdb.models['Model-1'].ViscoStep(cetol=1e-02, description='ViscoStep', 
				initialInc=maxincloading, maxInc=maxincloading, minInc=1e-5*maxincloading, maxNumInc=10*MinSteps, name='Loading', nlgeom=ON, previous='Initial', timePeriod=t_loading)
	elif Solvetype == 'Explicit':
		mdb.models['Model-1'].ExplicitDynamicsStep(improvedDtMethod=ON, 
			linearBulkViscosity=linbulk, maxIncrement=Time_exp/MinSteps, name='Loading', previous=
			'Initial', quadBulkViscosity=quadbulk, timePeriod=Time_exp)
		mdb.models['Model-1'].steps['Loading'].setValues(improvedDtMethod=ON, 
			massScaling=((DISABLE_THROUGHOUT_STEP, MODEL, THROUGHOUT_STEP, 0.0, 0.0, 
			None, 1, 0, 0.0, 0.0, 0, None), ))

	if UseStabilisation:
		mdb.models['Model-1'].steps['Loading'].setValues(adaptiveDampingRatio=0.05, 
			cetol=1e-05, continueDampingFactors=False, stabilizationMagnitude=2e-04, 
			stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
			

	del mdb.models['Model-1'].fieldOutputRequests['F-Output-1']
	del mdb.models['Model-1'].historyOutputRequests['H-Output-1']
	mdb.models['Model-1'].FieldOutputRequest(createStepName='Loading', name=
		'F-Output-1', timeInterval=t_loading/MinSteps , timeMarks=OFF, variables=('S', 'MISES', 
		'MISESMAX', 'TSHR', 'CTSHR', 'ALPHA', 'TRIAX', 'VS', 'PS', 'CS11', 
		'ALPHAN', 'SSAVG', 'MISESONLY', 'PRESSONLY', 'SEQUT', 'E', 'VE', 'PE', 
		'VEEQ', 'PEEQ', 'PEEQT', 'PEEQMAX', 'PEMAG', 'PEQC', 'EE', 'IE', 'CE', 
		'CEP', 'CEEQ', 'CEMAG', 'CESW', 'THE', 'NE', 'LE', 'TE', 'TEEQ', 'TEVOL', 
		'EEQUT', 'ER', 'SE', 'SPE', 'SEPE', 'SEE', 'SEP', 'SALPHA', 'U', 'UT', 
		'UR', 'V', 'VT', 'VR', 'RBANG', 'RBROT', 'RF', 'RT', 'RM', 'CF', 'SF', 
		'TF', 'VF', 'ESF1', 'NFORC', 'NFORCSO', 'RBFOR', 'BF', 'CORIOMAG', 
		'ROTAMAG', 'CENTMAG', 'CENTRIFMAG', 'GRAV', 'P', 'HP', 'TRSHR', 'TRNOR', 
		'ENER', 'ELEN', 'ELEDEN'))
	mdb.models['Model-1'].HistoryOutputRequest(createStepName='Loading', name=
		'H-Output-1', timeInterval=t_loading/MinSteps, variables=PRESELECT)

						
			
	if Dorelaxation:
		mdb.models['Model-1'].ViscoStep(cetol=1e-05, initialInc=0.01, maxInc=0.5*t_relaxend, 
			minInc=0.0001, maxNumInc=10000, name='Relaxation', previous='Loading', timePeriod=t_relaxend)
	
	
	#####################################
	# Constraints and Boundary Conditions if NOT rotated
	####################################


	# Reference Points
	#mdb.models['Model-1'].rootAssembly.ReferencePoint(point=((Nx-1)/2*d_cell, (Ny-0.)*d_cell, 0.0))
	#RPID1 = mdb.models['Model-1'].rootAssembly.referencePoints.keys()[0]

	if not  RotatedOrientation:



		if BCBottom =='support':
			mdb.models['Model-1'].rootAssembly.ReferencePoint(point=((Nx-1)/2*d_cell, -0.5*d_cell, 0.0))
			RPID2 = mdb.models['Model-1'].rootAssembly.referencePoints.keys()[0]





		#############
		RP_NxB = Nx*[None]


		for indX in range(Nx):
			
			mdb.models['Model-1'].rootAssembly.Set(name='BottomPoint-'+str(indX), vertices=
				mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].vertices.findAt(
				((d_cell*indX, 		-dh3, 0.0), ), ))

			mdb.models['Model-1'].rootAssembly.Set(edges=
				mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(
				((d_cell*indX-1e-5, 		-dh3, 0.0), ), ((d_cell*indX+1e-5, 		-dh3, 0.0), ), ), name='BottomEdges-'+str(indX))
			if UseConstraints :	
				mdb.models['Model-1'].Coupling(controlPoint=
					mdb.models['Model-1'].rootAssembly.sets['BottomPoint-'+str(indX)], couplingType=
					KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
					'BottomConstraint-'+str(indX), surface=mdb.models['Model-1'].rootAssembly.sets['BottomEdges-'+str(indX)], 
					u1=ON, u2=ON, ur3=OFF)	

		for indX in range(Nx):
			
			mdb.models['Model-1'].rootAssembly.Set(name='TopPoint-'+str(indX), vertices=
				mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].vertices.findAt(
				((d_cell*indX, 		(Ny-1)*d_cell+dh3, 0.0), ), ))

			mdb.models['Model-1'].rootAssembly.Set(edges=
				mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(
				((d_cell*indX-1e-5, 		(Ny-1)*d_cell+dh3, 0.0), ), ((d_cell*indX+1e-5, 		(Ny-1)*d_cell+dh3, 0.0), ), ), name='TopEdges-'+str(indX))
				
			if UseConstraints:	
				mdb.models['Model-1'].Coupling(controlPoint=
					mdb.models['Model-1'].rootAssembly.sets['TopPoint-'+str(indX)], couplingType=
					KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
					'TopConstraint-'+str(indX), surface=mdb.models['Model-1'].rootAssembly.sets['TopEdges-'+str(indX)], 
					u1=ON, u2=ON, ur3=OFF)		



		if BCBottom=='rigid':
			# Rigid side with constraint
			
			mdb.models['Model-1'].rootAssembly.Set(name='FixXPoint', vertices=
				mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].vertices.findAt(
				((-dh3, 0., 0.0), ), ))

			mdb.models['Model-1'].rootAssembly.Set(edges=
				mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(
				((-dh3, 1e-5, 0.0), ), ((-dh3, -1e-5, 0.0), ), ), name=		'FixXEdges')
				
			if UseConstraints:	
				mdb.models['Model-1'].Coupling(controlPoint=
					mdb.models['Model-1'].rootAssembly.sets['FixXPoint'], couplingType=
					KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
					'SideNodeConstraint', surface=mdb.models['Model-1'].rootAssembly.sets['FixXEdges'], 
					u1=ON, u2=ON, ur3=OFF)						
						

			if BCBottom =='support':
				mdb.models['Model-1'].Coupling(controlPoint=Region(referencePoints=(
					mdb.models['Model-1'].rootAssembly.referencePoints[RPID2], )), couplingType=
					KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
					'CouplingBottom', surface=
					mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].sets['Bottom']
					, u1=OFF, u2=ON, ur3=OFF)	

				mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
					distributionType=UNIFORM, fieldName='', localCsys=None, name='FixYBottom', 
					region=Region(referencePoints=(
					mdb.models['Model-1'].rootAssembly.referencePoints[RPID2], )), u1=SET, u2=SET, ur3=SET)
			elif BCBottom =='rigid':	
				mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
					distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
					'FixBottom', region=mdb.models['Model-1'].rootAssembly.instances['Rigid-2'].sets['RP_Rigid'], u1=SET, u2=
					SET, ur3=UNSET)
			if FixX:
				mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
					distributionType=UNIFORM, fieldName='', localCsys=None, name=
					'FixXBottom', region=Region(
					vertices=mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].vertices.findAt(
					((-dh3, 0., 0.0), ), )), u1=SET, u2=UNSET, ur3=UNSET)	
				
			if amptype == 'Straight':
				mdb.models['Model-1'].TabularAmplitude(data=((0.0, 0.0), (t_loading, 1.0)), name=
					'Amp-1', smooth=SOLVER_DEFAULT, timeSpan=TOTAL)
			elif amptype =='Smooth': 
				mdb.models['Model-1'].SmoothStepAmplitude(data=((0.0, 0.0), (t_loading, 1.0)), name=
					'Amp-1', timeSpan=STEP)	
				
				
			if UseLoadY:
				mdb.models['Model-1'].ConcentratedForce(amplitude='Amp-1', cf2=LoadY, 
					createStepName='Loading', distributionType=UNIFORM, field='', localCsys=
					None, name='Load-1', region=	mdb.models['Model-1'].rootAssembly.instances['Rigid-1'].sets['RP_Rigid']
					)
				mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Loading', 
					distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
					'FixXRigidTop', region=mdb.models['Model-1'].rootAssembly.instances['Rigid-1'].sets['RP_Rigid'], u1=0.0, u2=
					UNSET, ur3=UNSET)
				

			else:
				mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', createStepName='Loading', 
					distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
					'Compress', region=mdb.models['Model-1'].rootAssembly.instances['Rigid-1'].sets['RP_Rigid'], u1=0.0, u2=
					-d_cell, ur3=UNSET)
			
		if ImperfGrav:	
			mdb.models['Model-1'].Gravity(comp1=GravLoad, comp2=0.1*GravLoad, createStepName='Loading', 
				distributionType=UNIFORM, field='', name='GravLoad')






		###########################
		# Contact Definition
		###########################		
		if UseContact:
			Surfs =[]
			for indX in range(Nx):
				TX = (indX)*(d_cell)
				TY = (Ny-1)*d_cell+dh3
				Surfs = Surfs+[mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(((TX-1e-5,TY, 0.0), ))]
				Surfs = Surfs+[mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(((TX+1e-5,TY, 0.0), ))]
			Surfs = tuple(Surfs)
			mdb.models['Model-1'].rootAssembly.Surface(name='TopSurfaces', side2Edges=Surfs	)	
			mdb.models['Model-1'].rootAssembly.Set(edges=Surfs,name='TopContactSet'	)	
			mdb.models['Model-1'].rootAssembly.Surface(name='TopRigidSurface', side2Edges=
				mdb.models['Model-1'].rootAssembly.instances['Rigid-1'].edges.findAt(((
				0., (Ny-1)*d_cell+dh3+h_rigid, 0.0), )))




			
			if BCBottom =='rigid':
				Surfs =[]
				for indX in range(Nx):
					TX = (indX)*(d_cell)
					TY = -dh3
					Surfs = Surfs+[mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(((TX-1e-5,TY, 0.0), ))]
					Surfs = Surfs+[mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].edges.findAt(((TX+1e-5,TY, 0.0), ))]
				Surfs = tuple(Surfs)
				mdb.models['Model-1'].rootAssembly.Surface(name='BottomSurfaces', side2Edges=Surfs	)	
				mdb.models['Model-1'].rootAssembly.Set(edges=Surfs,name='BottomContactSet'	)	
				mdb.models['Model-1'].rootAssembly.Surface(name='BottomRigidSurface', side1Edges=
					mdb.models['Model-1'].rootAssembly.instances['Rigid-2'].edges.findAt(((
					0., -dh3-h_rigid, 0.0), )))		
				

			mdb.models['Model-1'].ContactProperty('ContactProp')
			mdb.models['Model-1'].interactionProperties['ContactProp'].TangentialBehavior(		formulation=FRICTIONLESS)	
			if FricCoef > 1e-10:
				mdb.models['Model-1'].interactionProperties['ContactProp'].tangentialBehavior.setValues(
					dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
					formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
					pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
					table=((FricCoef, ), ), temperatureDependency=OFF)		
			mdb.models['Model-1'].interactionProperties['ContactProp'].NormalBehavior(
				constraintEnforcementMethod=DEFAULT, contactStiffness=K_Contact, 
				pressureOverclosure=LINEAR)
			mdb.models['Model-1'].interactionProperties['ContactProp'].Damping(
				clearanceDependence=STEP, definition=DAMPING_COEFFICIENT, table=((alpha_rigid, ), 
				), tangentFraction=DEFAULT)		
				
			if Solvetype =='Implicit':	
				mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
					clearanceRegion=None, createStepName='Initial', datumAxis=None, 
					initialClearance=OMIT, interactionProperty='ContactProp', master=
					mdb.models['Model-1'].rootAssembly.surfaces['TopRigidSurface'], name=
					'SurfContact', slave=
					mdb.models['Model-1'].rootAssembly.sets['TopContactSet'], sliding=FINITE, 
					thickness=ON)
					
				if BCBottom=='rigid':
					mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
						clearanceRegion=None, createStepName='Initial', datumAxis=None, 
						initialClearance=OMIT, interactionProperty='ContactProp', master=
						mdb.models['Model-1'].rootAssembly.surfaces['BottomRigidSurface'], name=
						'SurfContactBot', slave=
						mdb.models['Model-1'].rootAssembly.sets['BottomContactSet'], sliding=FINITE, 
						thickness=ON)		
			elif Solvetype =='Explicit':
				mdb.models['Model-1'].SurfaceToSurfaceContactExp(clearanceRegion=None, createStepName='Initial', datumAxis=None, 
					initialClearance=OMIT, interactionProperty='ContactProp', master=
					mdb.models['Model-1'].rootAssembly.surfaces['TopRigidSurface'], name=
					'SurfContact', slave=
					mdb.models['Model-1'].rootAssembly.sets['TopContactSet'], sliding=FINITE)	
				if BCBottom=='rigid':
					mdb.models['Model-1'].SurfaceToSurfaceContactExp(clearanceRegion=None, createStepName='Initial', datumAxis=None, 
						initialClearance=OMIT, interactionProperty='ContactProp', master=
						mdb.models['Model-1'].rootAssembly.surfaces['BottomRigidSurface'], name=
						'SurfContactBot', slave=
						mdb.models['Model-1'].rootAssembly.sets['BottomContactSet'], sliding=FINITE)		


	#####################################
	# Constraints and Boundary Conditions if  rotated
	####################################


	if   RotatedOrientation:
		####################
		#Constraints

		RPTop = {}
		RPBottom={}
		if BCBottom =='support' and BCTop =='support':
		
			for indX in range(NBCrot):
				#Bottom side
				x=(DeltaBCrot+2.*indX)*a_cell
				y=-0.5*(Nyrot-1)*a_cell
				RPBottom[indX] = 	mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(x, y, 0.0))
				mdb.models['Model-1'].rootAssembly.Set(name='RPBottom-'+str(indX), referencePoints=(    mdb.models['Model-1'].rootAssembly.referencePoints[RPBottom[indX].id], ))
				mdb.models['Model-1'].rootAssembly.Set(faces=    mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].faces.findAt(    ((x, y, 0.0), )), name='FaceBottom-'+str(indX))
				mdb.models['Model-1'].Coupling(controlPoint=    mdb.models['Model-1'].rootAssembly.sets['RPBottom-'+str(indX)], couplingType=KINEMATIC,     influenceRadius=WHOLE_SURFACE, localCsys=None, name='CouplingBottom-'+str(indX),     surface=mdb.models['Model-1'].rootAssembly.sets['FaceBottom-'+str(indX)], u1=ON, u2=ON,     ur3=ON)
			for indX in range(NBCrot):
				#Top side
				x=(DeltaBCrot+2.*indX)*a_cell
				y=0.5*(Nyrot-1)*a_cell
				RPTop[indX] = 	mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(x, y, 0.0))
				mdb.models['Model-1'].rootAssembly.Set(name='RPTop-'+str(indX), referencePoints=(    mdb.models['Model-1'].rootAssembly.referencePoints[RPTop[indX].id], ))
				mdb.models['Model-1'].rootAssembly.Set(faces=    mdb.models['Model-1'].rootAssembly.instances['AssemblyPart-1'].faces.findAt(    ((x, y, 0.0), )), name='FaceTop-'+str(indX))
				mdb.models['Model-1'].Coupling(controlPoint=    mdb.models['Model-1'].rootAssembly.sets['RPTop-'+str(indX)], couplingType=KINEMATIC,     influenceRadius=WHOLE_SURFACE, localCsys=None, name='CouplingTop-'+str(indX),     surface=mdb.models['Model-1'].rootAssembly.sets['FaceTop-'+str(indX)], u1=ON, u2=ON,     ur3=ON)
			####################
			# Boundary conditions
			
			if amptype == 'Straight':
				mdb.models['Model-1'].TabularAmplitude(data=((0.0, 0.0), (t_loading, 1.0)), name=
					'Amp-1', smooth=SOLVER_DEFAULT, timeSpan=TOTAL)
			elif amptype =='Smooth': 
				mdb.models['Model-1'].SmoothStepAmplitude(data=((0.0, 0.0), (t_loading, 1.0)), name=
					'Amp-1', timeSpan=STEP)

			int_Fix = int(np.round(0.5*DeltaBCrot))-1
			for indX in range(NBCrot):
				if indX==int_Fix:
					mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=
						PERTURBATION_AND_BUCKLING, createStepName='Loading', distributionType=
						UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Fix_bottom', region=
						mdb.models['Model-1'].rootAssembly.sets['RPBottom-'+str(indX)], u1=0.0, u2=0.0, ur3=
						UNSET)
					if LoadingWay == 'Force'	:
					

						mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=
							PERTURBATION_AND_BUCKLING, createStepName='Loading', distributionType=
							UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Fix_x_top', region=
							mdb.models['Model-1'].rootAssembly.sets['RPTop-'+str(indX)], u1=0.0, u2=UNSET, ur3=
							UNSET)

				else:
					if LoadingWay == 'Force':
						mdb.models['Model-1'].ConcentratedForce(cf2=-LoadY/NBCrot*(1.+Imperfload*indX/NBCrot), createStepName='Loading', 
							distributionType=UNIFORM, field='', localCsys=None, name='Compress_bottom-'+str(indX), 
							region=mdb.models['Model-1'].rootAssembly.sets['RPBottom-'+str(indX)])	
			
					elif LoadingWay =='DispImp':
						
						mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', buckleCase=
							PERTURBATION_AND_BUCKLING, createStepName='Loading', distributionType=
							UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Fix_y_bottom-'+str(indX), region=
							mdb.models['Model-1'].rootAssembly.sets['RPBottom-'+str(indX)], u1=UNSET, u2=0.0, ur3=
							UNSET)
							
				if LoadingWay=='Force':
					mdb.models['Model-1'].ConcentratedForce(cf2=LoadY/NBCrot*(1.+Imperfload*indX/NBCrot), createStepName='Loading', 
						distributionType=UNIFORM, field='', localCsys=None, name='Compress_top-'+str(indX), 
						region=mdb.models['Model-1'].rootAssembly.sets['RPTop-'+str(indX)])
				
				elif LoadingWay =='DispImp':
					mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', createStepName='Loading', 
						distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
						'Compress_top-'+str(indX), region=mdb.models['Model-1'].rootAssembly.sets['RPTop-'+str(indX)], u1=UNSET, u2=
						DispY*(1.+Imperfload*indX/NBCrot), ur3=UNSET)			


	###########################
	# Job
	###########################		
	
	mdb.models['Model-1'].rootAssembly.regenerate()

	
		
	if Dorelaxation:	
		mdb.models['Model-1'].TimePoint(name='TimePoints-1', points=t_relaxtuple )
		mdb.models['Model-1'].FieldOutputRequest(createStepName='Relaxation', name=		'F-Output-2', timePoint='TimePoints-1', variables=('S', 'MISES', 
			'MISESMAX', 'TSHR', 'CTSHR', 'ALPHA', 'TRIAX', 'VS', 'PS', 'CS11', 
			'ALPHAN', 'SSAVG', 'MISESONLY', 'PRESSONLY', 'SEQUT', 'E', 'VE', 'PE', 
			'VEEQ', 'PEEQ', 'PEEQT', 'PEEQMAX', 'PEMAG', 'PEQC', 'EE', 'IE', 'CE', 
			'CEP', 'CEEQ', 'CEMAG', 'CESW', 'THE', 'NE', 'LE', 'TE', 'TEEQ', 'TEVOL', 
			'EEQUT', 'ER', 'SE', 'SPE', 'SEPE', 'SEE', 'SEP', 'SALPHA', 'U', 'UT', 
			'UR', 'V', 'VT', 'VR', 'RBANG', 'RBROT', 'RF', 'RT', 'RM', 'CF', 'SF', 
			'TF', 'VF', 'ESF1', 'NFORC', 'NFORCSO', 'RBFOR', 'BF', 'CORIOMAG', 
			'ROTAMAG', 'CENTMAG', 'CENTRIFMAG', 'GRAV', 'P', 'HP', 'TRSHR', 'TRNOR', 
			'ENER', 'ELEN', 'ELEDEN'))
		mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].deactivate(
			'Relaxation')
			
		mdb.models['Model-1'].HistoryOutputRequest(createStepName='Relaxation', name='H-Output-2', timeMarks=OFF, timePoint='TimePoints-1', variables=PRESELECT)
		mdb.models['Model-1'].historyOutputRequests['H-Output-1'].deactivate(
			'Relaxation')		
	""" 
	mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
		'S', 'MISES', 'MISESMAX', 'TSHR', 'CTSHR', 'ALPHA', 'TRIAX', 'VS', 'PS', 
		'CS11', 'SSAVG', 'MISESONLY', 'PRESSONLY', 'SEQUT', 'E', 'VE', 'PE', 
		'VEEQ', 'PEEQ', 'PEEQT', 'PEEQMAX', 'PEMAG', 'PEQC', 'EE', 'IE', 'THE', 
		'NE', 'LE', 'TE', 'TEEQ', 'TEVOL', 'EEQUT', 'ER', 'SE', 'SPE', 'SEPE', 
		'SEE', 'SEP', 'SALPHA', 'U', 'RF', 'RT', 'RM', 'CF', 'SF', 'TF', 'VF', 
		'ESF1', 'NFORC', 'NFORCSO', 'RBFOR'))		
	""" 


	#compressive buckling
	jobtitle = jobtitlestart+'-Mat-'+str(ind_ItMat)
	myjob = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
		explicitPrecision=PrecisionLevelExp, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
		multiprocessingMode=DEFAULT, name=jobtitle			, nodalOutputPrecision=SINGLE, 
		numCpus=4, numGPUs=0, numDomains=4, queue=None, resultsFormat=ODB, scratch='', type=
		ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)	

	


	
	if SubmitJob or writeInput:
		time.sleep(3)
		if writeInput:
			myjob.writeInput(		consistencyChecking=OFF)

		
		if SubmitJob:
			try:
				myjob.submit(consistencyChecking=	OFF)
				myjob.waitForCompletion()	
			except:
				pass
			time.sleep(3)
		
		
		if N_ItMat>1.1:
			for ind_ItMat in range(1,N_ItMat):
				Mat  = Materials[2]				
				E_TriH = E_TriH_vec[ind_ItMat]
				lambda_TriH = E_TriH*nu_H/((1+nu_H)*(1-2*nu_H))
				mu_TriH=E_TriH/(2*(1+nu_H))
				C10_TriH = mu_TriH/2
				D1_TriH = lambda_TriH/2
				
				
				if MatPot =='Linear':
					mdb.models['Model-1'].materials[Mat].elastic.setValues(table=((E_TriH, nu_H), ))
				elif MatPot == 'NeoHooke':	
					mdb.models['Model-1'].materials[Mat].hyperelastic.setValues(table=((C10_TriH, D1_TriH), ))
				elif MatPot == 'ArrudaBoyce':	
					mdb.models['Model-1'].materials[Mat].hyperelastic.setValues(table=((2*C10_TriH, lam_AB, D1_TriH), ))
				elif MatPot == 'ArrudaBoyceVisco':	
					mdb.models['Model-1'].materials[Mat].hyperelastic.setValues(table=((2*C10_TriH, lam_AB, D1_TriH), ))
					Mat  = Materials[1]				
					mdb.models['Model-1'].materials[Mat].viscoelastic.setValues(
						domain=TIME, table=Pronytermsvisco[ind_ItMat], time=PRONY)
				
				for Mat in Materials:
					mdb.models['Model-1'].materials[Mat].density.setValues(table=((density_vec[ind_ItMat], ), ))
				#Loading				
				t_loading = t_loading_vec[ind_ItMat]
				
				if steptype =='ViscoStep':
					t_loading = t_loading_vec[ind_ItMat]
					maxincloading = t_loading/MinSteps 
					mdb.models['Model-1'].steps['Loading'].setValues(cetol=1e-02, 
						initialInc=maxincloading, maxInc=maxincloading, minInc=1e-3*maxincloading, maxNumInc=10*MinSteps, timePeriod=t_loading)
					mdb.models['Model-1'].amplitudes['Amp-1'].setValues(data=((0.0, 0.0), (t_loading, 
						1.0)), timeSpan=STEP)


					del mdb.models['Model-1'].fieldOutputRequests['F-Output-1']
					del mdb.models['Model-1'].historyOutputRequests['H-Output-1']
					mdb.models['Model-1'].FieldOutputRequest(createStepName='Loading', name=
						'F-Output-1', timeInterval=t_loading/MinSteps , timeMarks=OFF, variables=('S', 'MISES', 
						'MISESMAX', 'TSHR', 'CTSHR', 'ALPHA', 'TRIAX', 'VS', 'PS', 'CS11', 
						'ALPHAN', 'SSAVG', 'MISESONLY', 'PRESSONLY', 'SEQUT', 'E', 'VE', 'PE', 
						'VEEQ', 'PEEQ', 'PEEQT', 'PEEQMAX', 'PEMAG', 'PEQC', 'EE', 'IE', 'CE', 
						'CEP', 'CEEQ', 'CEMAG', 'CESW', 'THE', 'NE', 'LE', 'TE', 'TEEQ', 'TEVOL', 
						'EEQUT', 'ER', 'SE', 'SPE', 'SEPE', 'SEE', 'SEP', 'SALPHA', 'U', 'UT', 
						'UR', 'V', 'VT', 'VR', 'RBANG', 'RBROT', 'RF', 'RT', 'RM', 'CF', 'SF', 
						'TF', 'VF', 'ESF1', 'NFORC', 'NFORCSO', 'RBFOR', 'BF', 'CORIOMAG', 
						'ROTAMAG', 'CENTMAG', 'CENTRIFMAG', 'GRAV', 'P', 'HP', 'TRSHR', 'TRNOR', 
						'ENER', 'ELEN', 'ELEDEN'))
					mdb.models['Model-1'].HistoryOutputRequest(createStepName='Loading', name=
						'H-Output-1', timeInterval=t_loading/MinSteps, variables=PRESELECT)

							
						
					if Dorelaxation:	
						del mdb.models['Model-1'].fieldOutputRequests['F-Output-2']
						del mdb.models['Model-1'].historyOutputRequests['H-Output-2']
						mdb.models['Model-1'].FieldOutputRequest(createStepName='Relaxation', name=		'F-Output-2', timePoint='TimePoints-1', variables=('S', 'MISES', 'MISESMAX', 'TSHR', 'CTSHR', 'ALPHA', 'TRIAX', 'VS', 'PS', 'CS11', 
							'ALPHAN', 'SSAVG', 'MISESONLY', 'PRESSONLY', 'SEQUT', 'E', 'VE', 'PE', 
							'VEEQ', 'PEEQ', 'PEEQT', 'PEEQMAX', 'PEMAG', 'PEQC', 'EE', 'IE', 'CE', 
							'CEP', 'CEEQ', 'CEMAG', 'CESW', 'THE', 'NE', 'LE', 'TE', 'TEEQ', 'TEVOL', 
							'EEQUT', 'ER', 'SE', 'SPE', 'SEPE', 'SEE', 'SEP', 'SALPHA', 'U', 'UT', 
							'UR', 'V', 'VT', 'VR', 'RBANG', 'RBROT', 'RF', 'RT', 'RM', 'CF', 'SF', 
							'TF', 'VF', 'ESF1', 'NFORC', 'NFORCSO', 'RBFOR', 'BF', 'CORIOMAG', 
							'ROTAMAG', 'CENTMAG', 'CENTRIFMAG', 'GRAV', 'P', 'HP', 'TRSHR', 'TRNOR', 
							'ENER', 'ELEN', 'ELEDEN'))
						mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].deactivate(
							'Relaxation')
							
						mdb.models['Model-1'].HistoryOutputRequest(createStepName='Relaxation', name='H-Output-2', timeMarks=OFF, timePoint='TimePoints-1', variables=PRESELECT)
						mdb.models['Model-1'].historyOutputRequests['H-Output-1'].deactivate(
							'Relaxation')		
						
						
				jobtitleold = jobtitle
				jobtitle = jobtitlestart+'-Mat-'+str(ind_ItMat)
				mdb.jobs.changeKey(fromName=jobtitleold, toName=
				jobtitle)
				myjob = mdb.jobs[jobtitle]
				#mdb.saveAs(pathName=file_name)	
				time.sleep(3)
						
				if writeInput:
					myjob.writeInput(		consistencyChecking=OFF)
				if SubmitJob:
					try:
						myjob.submit(consistencyChecking=	OFF)
						myjob.waitForCompletion()	
					except:
						pass
				time.sleep(3)
	
				

	mdb.saveAs(pathName=file_name)	








#################
#################
#################
#################
#################
#################
#################
#################
#################
#################
#################
#################
#################



# FOR NEXT TIME: ADD CONTACT OR RIGID ELEMENTS ON TOP AND BOTTOM!!!