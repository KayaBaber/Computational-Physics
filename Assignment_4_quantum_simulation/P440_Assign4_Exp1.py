'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 4 - Quantum Simulations
Exploration 1
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt



#make x array of N intervals, without first and last points, so N-1 points
#make banded matrix using diag and 2nd order central difference coeficients repeated, don't do edges
#take eiganenergies and eiganvectors from banded matrix
#add first and last values to x and eiganvectors (0)
#plot first 8 eiganmodes and some spurious nonphysical modes
#look at 8th eiganenergy for different N


def make_banded(N):
    N-=2
    stepSize = 1.0/N
    bandTopBot = [1.0/ (stepSize**2)]*(N-1)
    bandMid = [-2.0/ (stepSize**2)]*N
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded

N=100
banded = make_banded(N)

xArray=np.linspace(-0.5,0.5,N)
xArraySub=xArray[1:-1]
#print xArraySub
print banded
eiganVal, eiganVect=LA.eig(banded)
print eiganVect
for i in range(25):
    psi=np.insert(eiganVect[:,i],0 ,0)
    psi=np.append(psi,0)
    #print eiganVal
    #print eiganVect[0]
    #print psi
    plt.plot(xArray, psi)
    plt.grid()
    plt.savefig('plots/eiganVect'+str(i)+'.png',bbox_inches='tight')
    #plt.show()
    plt.clf()


#energyArray=[]
#for i in range(15,50,1):
#    banded = make_banded(i)
#    energyArray.append(LA.eig(banded)[0][8])
#plt.plot(range(15,50,1),energyArray)
#plt.show()
    
    
