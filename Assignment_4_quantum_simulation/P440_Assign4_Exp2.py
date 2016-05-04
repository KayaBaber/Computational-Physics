'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 4 - Quantum Simulations
Exploration 1
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

#continue from exploration 1
#normalize eigan vectors
    #square eigan vector points
    #use trap rule to integrate over x
        #sum up all points
        #divide by 2 times number of intervals 2(N-1)
    #take square root of integral
    #divide all eigan vector points by the square root
#integrate over x to get delta_E
    #square normalized eigan vector points
    #multiple squared points by V'(x)
        #V'(x) = V_0*exp(-256*x^2)
        #V_0 = E_0 / 10.0
    #use trap rule to integrate over x
        #sum up all points
        #divide by 2 times number of intervals 2(N-1)
#print out the original eigan energies
#print out the delta_E values
#add the delta_E values to the original eigan energies
#print out the new eigan energies


def make_banded(N):
    N-=2    #we reduce N by 2 because we aren't using the outer ring of indicies
    stepSize = 1.0/N
    bandTopBot = [1.0/ (stepSize**2)]*(N-1)
    bandMid = [-2.0/ (stepSize**2)]*N
    
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded


def normalize(psi):
    psiSqr = psi**2
    psiSqrSum = np.sum(psiSqr)
    normConst = psiSqrSum ** (1./2.)
    psiNorm = psi / normConst
    return psi
    
    

N=100
banded = make_banded(N)

xArray=np.linspace(-0.5,0.5,N)
xArraySub=xArray[1:-1]

eiganVal, eiganVect=LA.eig(banded)

num_vectors = 25

psiArray=[]
deltaEArray=[]
for i in range(num_vectors):
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    psi=normalize(psi)
    psiArray.append(psi)
    


