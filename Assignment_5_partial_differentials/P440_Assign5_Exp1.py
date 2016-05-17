'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 5 - PDEs
Exploration 1 - Parabolic PDEs: Thermal Diffusion
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
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


N = 10          #number of spatial steps per skin depth       
M = 10         #number of temporal steps per period
periods = 2   #number of periods
skins = 10      #number of skin depths


operatorCrank = make_operator(N*skins,M)

thermal = [1.]*N*skins

#simulation loop
thermLog = []
for i in range(M*periods):
    plt.plot(np.linspace(0,skins,N*skins),thermal)
    plt.show()    
    thermLog.append(thermal)
    thermal = np.dot(operatorCrank,thermal)
    thermal[0] = 1. + math.sin(2.*math.pi*i)
    thermal[-1] = 1.

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.linspace(0,skins,N*skins)
y = np.linspace(0,periods,M*periods)
XX, YY = np.meshgrid(x,y)
ZZ = thermLog

ax.plot_surface(XX, YY, ZZ, 
                #cmap=cm.spectral
               )
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Temperature')
plt.show()
