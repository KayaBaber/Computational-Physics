'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 5 - PDEs
Exploration 2 - Parabolic PDEs: The Wave Equation
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math
import cmath



L = 2.*math.pi  #set the x range to (0->2pi)
N = 1000        #number of spatial intervals and points (since it loops)
steps = 1000    #number of timesteps
stepSize = 0.1  #temporal step size

#make initial velocity array in real space
velPhysGauss = np.exp(-10 * ( (np.linspace(0,L,N+1)[:-1]-math.pi) ** 2)) + [0J]*N

#make initial density array in real space
denPhysFlat = [0 + 0J]*N
print denPhysFlat
print denPhysFlat[2]+2.
#fft both to fourier space
#make a column vector of density_f appended to velocity_f

#make derivative operator matrix
    #i(diag (0,1,2,3,4,.. -4, -3, -2, -1)), where i is imaginary i
#make quad matrix [[I][op],[op][I]] and negative quad [[I][-op],[-op][I]]
#step forward
    #matrix multiply the negative quad by the FFT column vector
    #linear algebra solve the pos_quad*newFFTvector = above-result for newFFTvector
    #make seperate copies of the velocity_f and density_f components of newFFTvector
    #inverse FFT the components and append the real parts of velocity and density to a log
#repeat stepping for num steps
#plot the velocity and density logs in 3D
#maybe make an animation





 
