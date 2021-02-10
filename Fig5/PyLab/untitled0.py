# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:59:55 2019

@author: David
"""

import numpy as np



a10=np.array([[0,1,2,3],
            [1,0,1,2],
            [2,1,0,1],
            [3,2,1,0]])


a1=np.zeros((4,4))
for ind in range(4):
    a1[ind,:]=a10[ind,:]*4/np.sum(a10[ind,:])
    
    
    
    
    



a2= np.ones((4,4))

b1=1/2*a1**2
b2=1/2*a2**2


c1=np.sum(b1)
c2=np.sum(b2)

    
