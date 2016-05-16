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

def make_quad_ops(N):
    iden = np.identity(N)
    #make derivative operator matrix
    #i(diag (0,1,2,3,4,.. -4, -3, -2, -1)), where i is imaginary i
    opFlat = np.append(np.linspace(0,N/2, (N/2)+1)*1J,
                       np.linspace((-N/2)+1,0,(N/2))[:-1]*1J)
    opDiag = np.diag(opFlat)
    top = np.append(iden, opDiag, axis=1)
    bot = np.append(opDiag, iden, axis=1)
    posQuad = np.append(top, bot, axis=0)
    negTop = np.append(iden, -opDiag, axis=1)
    negBot = np.append(-opDiag, iden, axis=1)
    negQuad = np.append(negTop, negBot, axis=0)
    return posQuad, negQuad


def step_forward(velDen,posQuad,negQuad):
    #matrix multiply the negative quad by the FFT column vector
    RHS = np.dot(negQuad,velDen)
    #linear algebra solve the pos_quad*newFFTvector = above-result for newFFTvector
    solution = LA.tensorsolve(posQuad,RHS)
    return solution



L = 2.*math.pi  #set the x range to (0->2pi)
N = 10        #number of spatial intervals and points (since it loops)
steps = 1000    #number of timesteps
stepSize = 0.1  #temporal step size

#make initial velocity array in physical space
velPhysGauss = np.exp(-10 * ( (np.linspace(0,L,N+1)[:-1]-math.pi) ** 2)) + [0J]*N

#make initial density array in physical space
denPhysFlat = [0 + 0J]*N

#fft both to fourier space
velF = np.fft.fft(velPhysGauss)
denF = np.fft.fft(denPhysFlat)

#make a column vector of density_f appended to velocity_f
velDen = velF
velDen = np.append(velDen,denF)

#make quad matrix [[I][op],[op][I]] and negative quad [[I][-op],[-op][I]]
posQuad, negQuad = make_quad_ops(N)

#step forward
velDen = step_forward(
#make seperate copies of the velocity_f and density_f components of newFFTvector
#inverse FFT the components and append the real parts of velocity and density to a log
#repeat stepping for num steps
#plot the velocity and density logs in 3D
#maybe make an animation





 
