'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 4 - Quantum Simulations
Exploration 1
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math

def make_banded(N):
    N-=2                                 #we reduce N by 2 because we aren't using the outer ring of indicies
    stepSize = 1.0/N  
    bandTopBot = [1.0/ (stepSize**2)]*(N-1)
    bandMid = [-2.0/ (stepSize**2)]*N
    
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded
#continue from exploration 1


#normalize eigan vectors
def normalize(psi):
    #takes in eigan vector and normalizes it
    psiSqr = psi**2                     #square eigan vector points
    #use trap rule to integrate over x
    psiSqrSum = np.sum(psiSqr) / (2.* (len(psi)-1) ) #sum up all points and divide by 2 times number of intervals 2(N-1)
    normConst = psiSqrSum ** (1./2.)    #take square root of integral
    psiNorm = psi / normConst           #divide all eigan vector points by the square root
    return psiNorm


def V(x,V0):                            #potential function
    #V'(x) = V_0*exp(-256*x^2)
    #V_0 = E_0 / 10.0
    vOut = V0*math.exp(-256*(x**2))
    return vOut


def perturb(psi,x,V,V0):
    #integrate over x to get delta_E
    psiSqr = psi**2                     #square normalized eigan vector points
    for i in range(len(psi)):                   
        psiSqr[i] = psiSqr[i]*V(x[i],V0)    #multiply squared points by V'(x)
    #use trap rule to integrate over x
    psiSqrSum = np.sum(psiSqr) / (2.* (len(psi)-1) ) #sum up all points and divide by 2 times number of intervals 2(N-1)
    deltaE=psiSqrSum                    #the change in energy is the integral
    return deltaE


N=100
banded = make_banded(N)
x=np.linspace(-0.5,0.5,N)               #positions
xArraySub=x[1:-1]

eiganVal, eiganVect=LA.eig(banded)      #compute eigan vectors and values

numVectors = 25

psiArray=[]
deltaEArray=[]
for i in range(numVectors):             #constructs array of eigan vectors
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    psi=normalize(psi)
    psiArray.append(psi)
    
psiNormArray=[]
for p in psiArray:
    psiNorm = normalize(p)
    psiNormArray.append(psiNorm)
    
    
    

V0=eiganVal[2] / 10.0
deltaEArray
for p in psiNormArray:
    deltaE=perturb(p,x,V,V0)
    deltaEArray.append(deltaE)
print eiganVal[:numVectors]
newEnergy = np.add(deltaEArray,eiganVal[:numVectors])
print np.array(newEnergy)


vArray=[]
for i in xArraySub:
    vArray.append(V(i,V0))
Vdiag=np.diag(vArray)
H = np.add(Vdiag,banded)
eiganValH, eiganVect=LA.eig(H)
psiArray=[]
for i in range(numVectors):             #constructs array of eigan vectors
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    psi=normalize(psi)
    psiArray.append(psi)
    plt.plot(x, psi)
    plt.ylabel("Eigan Vector")
    plt.xlabel("Position (x)")
    plt.grid()
    plt.savefig('Assignment_4_quantum_simulation/plots/2AeiganVect'+str(i)+'.png',bbox_inches='tight')
    #plt.show()
    plt.clf()




V0=eiganVal[2]
deltaEArray=[]
for p in psiNormArray:
    deltaE=perturb(p,x,V,V0)
    deltaEArray.append(deltaE)
newEnergy = np.add(deltaEArray,eiganVal[:numVectors])
print np.array(newEnergy)

vArray=[]
for i in xArraySub:
    vArray.append(V(i,V0))
Vdiag=np.diag(vArray)
H = np.add(Vdiag,banded)
eiganValH, eiganVect=LA.eig(H)
psiArray=[]
for i in range(numVectors):             #constructs array of eigan vectors
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    psi=normalize(psi)
    psiArray.append(psi)
    plt.plot(x, psi)
    plt.ylabel("Eigan Vector")
    plt.xlabel("Position (x)")
    plt.grid()
    plt.savefig('Assignment_4_quantum_simulation/plots/2BeiganVect'+str(i)+'.png',bbox_inches='tight')
    #plt.show()
    plt.clf()
