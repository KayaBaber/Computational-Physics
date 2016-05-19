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
    return bandedCrank, negativeCrank
    

def step_forward(thermal,posCrank,negCrank):
    #matrix multiply the negative quad by the FFT column vector
    RHS = np.dot(negCrank,thermal)
    #linear algebra solve the pos_quad*newFFTvector = above-result for newFFTvector
    newThermal = LA.tensorsolve(posCrank, RHS)
    return newThermal
    
    
N = 100          #number of spatial steps per skin depth       
M = 100         #number of temporal steps per period
periods = 2     #number of periods
skins = 10      #number of skin depths


posCrank, negCrank = make_operator(N*skins,M)

thermal = [1.]*N*skins

#simulation loop
topLog = []
thermLog = []
for i in range(M*periods):
#    plt.plot(np.linspace(0,skins,N*skins),thermal)
#    plt.show()    
    thermLog.append(thermal)
    thermal = step_forward(thermal,posCrank,negCrank)
    thermal[0] = 1.0 + math.sin(2.0*math.pi*(i/float(M)))
    topLog.append(thermal[0])
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

#plt.plot(np.linspace(0,periods,M*periods),topLog)
#plt.show() 