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
#make a negative copy of the matrix
#replace first and last rows of the first matrix with [1,0,0,0,0...0] and [0,0,0,0,0,...0,1]
#invert first matrix
#multiply the two matrices
#initialize column vector 
#matrix multiply the final matrix with the column vector
#replace boundaries of the column vector with [1+sin(2*pi*t)] and [1]
#repeat for number of timesteps


def make_banded(N,M):
    bandTopBot = [-1.]*(N-1)
    bandMid = [2. + (4.*M)/(N**2) ]*N
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded
    
    
def make_operator(N,M):
    bandedCrank = make_banded(N,M)
    negativeCrank = bandedCrank * (-1)
    bandedCrank[0] = [1] + [0]*(N-1)
    bandedCrank[-1] = [0]*(N-1) + [1]
    invertedCrank = LA.inv(bandedCrank)
    operatorCrank = invertedCrank.dot(negativeCrank)
    return operatorCrank

def operator(opertorCrank,thermal):
    
    
    return newThermal

N = 10          #number of spatial steps per skin depth       
M = 100         #number of temporal steps per period
periods = 100   #number of periods
skins = 10      #number of skin depths


operatorCrank = make_operator(N,M)
thermal = [1.]*N*skins
print thermal
#banded.dot(thermalArray)
