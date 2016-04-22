'''
Kaya Baber
Physics 440 - Computational Physics
Assignment 4 - Quantum Simulations
Exploration 1
'''

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


#make time array of N intervals, without first and last points
#make banded matrix using diag and 2nd order central difference coeficients repeated, don't do edges
#take eiganenergies and eiganvectors from banded matrix
#add first and last values to time and eiganvectors (0)
#plot first 8 eiganmodes and some spurious nonphysical modes
#look at 8th eiganenergy for different N


def make_banded(N):
    N-=2
    bandTopBot = [1]*(N-1)
    bandMid = [-2]*N
    banded = np.diag(bandMid)
    banded = np.add(banded,np.diag(bandTopBot,1))
    banded = np.add(banded,np.diag(bandTopBot,-1))
    return banded
    
banded = make_banded(25)
print banded
print LA.eig(banded)[0]    
energyArray=[]
for i in range(15,215,10):
    banded = make_banded(i)
    energyArray.append(LA.eig(banded)[0][8])
plt.plot(range(15,215,10),energyArray)
plt.show()
    
    
