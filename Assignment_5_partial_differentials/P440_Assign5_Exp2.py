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

#make initial velocity array in real space
#make initial density array in real space
#fft both to fourier space
#make a column vector of density_f appended to velocity_f
#set the x range to (0->2pi)
#make derivative operator matrix
    #i(diag (0,1,2,3,4,.. -4, -3, -2, -1)), where i is imaginary i
#make quad matrix [[I][op],[op][I]] and negative quad [[I][-op],[-op][I]]
#matrix multiply the negative quad by the FFT column vector
#linear algebra solve the pos_quad*newFFTvector = above-result for newFFTvector
#make seperate copies of the 
#inverse FFT newFFTvector and append the real parts of velocity and density
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


 
