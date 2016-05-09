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
    stepSize = 1.0/N 
    N-=2                                 #we reduce N by 2 because we aren't using the outer ring of indicies
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


N=65
banded = make_banded(N)
x=np.linspace(-0.5,0.5,N)               #positions
xArraySub=x[1:-1]

eiganVal, eiganVect=LA.eig(banded)      #compute eigan vectors and values

numVectors = 50

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


for i in range(40,numVectors):
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    plt.plot(x, psi)
    plt.ylabel("Psi")
    plt.xlabel("Position (x)")
    plt.title("Unperturbed")
    plt.grid()
    plt.savefig('Assignment_4_quantum_simulation/plots/1eiganVect'+str(i)+'.png',bbox_inches='tight')
    #plt.show()
    plt.clf()
    
 
V0=eiganVal[40] / 10.0          #GROUND STATE CHANGES WITH N, RECHECK EACH TIME
deltaEArray
for p in psiNormArray:
    deltaE=perturb(p,x,V,V0)
    deltaEArray.append(deltaE)
newEnergy = np.add(deltaEArray,eiganVal[:numVectors])



vArray=[]
for i in xArraySub:
    vArray.append(V(i,V0))
Vdiag=np.diag(vArray)
H = np.add(Vdiag,banded)
eiganValH, eiganVect=LA.eig(H)
eiganVect[:,40] = - eiganVect[:,40]
psiArray=[]
for i in range(40,numVectors):             #constructs array of eigan vectors
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    psi=normalize(psi)
    psiArray.append(psi)
    plt.plot(x, psi)
    plt.ylabel("Psi")
    plt.xlabel("Position (x)")
    plt.title("V_0 = 10% of Ground State Energy")
    plt.grid()
    plt.savefig('Assignment_4_quantum_simulation/plots/2AeiganVect'+str(i)+'.png',bbox_inches='tight')
    #plt.show()
    plt.clf()




V0=eiganVal[40]                 #GROUND STATE CHANGES WITH N, RECHECK EACH TIME
deltaEArray=[]
for p in psiNormArray:
    deltaE=perturb(p,x,V,V0)
    deltaEArray.append(deltaE)
newEnergy2 = np.add(deltaEArray,eiganVal[:numVectors])


vArray=[]
for i in xArraySub:
    vArray.append(V(i,V0))
Vdiag=np.diag(vArray)
H = np.add(Vdiag,banded)
eiganValH2, eiganVect=LA.eig(H)
psiArray=[]
for i in range(40,numVectors):             #constructs array of eigan vectors
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    psi=normalize(psi)
    psiArray.append(psi)
    plt.plot(x, psi)
    plt.ylabel("Psi")
    plt.xlabel("Position (x)")
    plt.title("V_0 = Ground State Energy")
    plt.grid()
    plt.savefig('Assignment_4_quantum_simulation/plots/2BeiganVect'+str(i)+'.png',bbox_inches='tight')
    #plt.show()
    plt.clf()
    
    
    
print "V0 = 10 % of Groundstate\nPertubation Theory:"
print np.array(newEnergy[40:])
print "Computed:"
print eiganValH[40:50]
print "Difference:"
print np.substract(np.array(newEnergy[40:]),eiganValH[40:50])
print "V0 = Groundstate\nPertubation Theory:"
print np.array(newEnergy2[40:])
print "Computed:"
print eiganValH2[40:50]
print "Difference:"
print np.subtract(np.array(newEnergy2[40:]),eiganValH2[40:50])
