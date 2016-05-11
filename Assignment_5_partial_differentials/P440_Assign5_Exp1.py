'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 5 - PDEs
Exploration 1 - Parabolic PDEs: Thermal Diffusion
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math




#make banded matrix 
#initialize column vector 
#matrix multiply
#modifiy boundaries
#repeat



def make_banded(N,M):
    bandTopBot = [2 + (4*M)/(N**2) ]*(N-1)
    bandMid = [-2.0/ (stepSize**2)]*N
    
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded
    
    
#banded.dot(thermalArray)
